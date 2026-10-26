#===============================================================================
# MS Thesis
#===============================================================================

# MUNICIPALITY level instead of the state level, extends the years for the DID, and adds a binary drought shock.

#===============================================================================
# 1). Preliminary 
#===============================================================================

rm(list=ls())
# Run time
start_time <- Sys.time()   

# PROJ fix
Sys.setenv(PROJ_LIB  = system.file("proj", package = "sf"),
           PROJ_DATA = system.file("proj", package = "sf"))

# Load packages
options(repos = c(CRAN = "https://cloud.r-project.org"))
packagelist <- c("terra","sf","data.table","readxl","haven")

newpackages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(newpackages)) install.packages(newpackages)

sapply(packagelist, function(i) require(i, character.only=TRUE))
rm(newpackages, packagelist)
# Progress bar
terraOptions(progress = 1)   


# Directories
dir <- list()
dir$root <- dirname(getwd()) 
# If the project root is not one folder up from the working folder, go up from the working folder until the folder that has Scripts/00_Master.do
if (!dir.exists(file.path(dir$root, "Data"))) {
  dir$root <- getwd()
  while (!file.exists(file.path(dir$root, "Scripts", "00_Master.do")) && dirname(dir$root) != dir$root) dir$root <- dirname(dir$root)
  if (!file.exists(file.path(dir$root, "Scripts", "00_Master.do"))) stop("Project folder not found. Set the working directory to the Scripts folder (open MS Thesis.Rproj) or run 00_Master.do")
}
dir$raw <- file.path(dir$root, "Data", "Raw")
dir$clean <- file.path(dir$root, "Data", "Cleaned")
dir$output_figure <- file.path(dir$root, "Output", "Figures")
# Create the folders this script writes to (does nothing if they already exist)
invisible(lapply(c(file.path(dir$clean, "Rainfall"), dir$output_figure), dir.create, recursive = TRUE, showWarnings = FALSE)) 

#===============================================================================
# Import rainfall data
#===============================================================================

# import the raster
chirps <- rast(file.path(dir$raw, "Rainfall", "chirps-v2.0.monthly.nc"))

# Create time index monthly and then subset the raster chirps to these dates, baseline is 1981-2004 which is 24 years 
time_index <- seq(as.Date("1981-01-01"), by="month", length.out=nlyr(chirps))
# We set the end to 2013 but techincally with DID shock it could be sooner but it doesn't hurt
keep <- which(time_index >= as.Date("1981-01-01") & time_index <= as.Date("2013-12-01"))
chirps <- chirps[[keep]]
# Rename layers from precip_121...precip_384 to their dates, e.g. "1981-01", so extracted columns are readable
names(chirps) <- format(time_index[keep], "%Y-%m")

# Import the shape file of the borders of mexico and then crop and mask will limit the raster to mexico
# INEGI uses a Lambert projection, so it is projected to the raster's lon/lat first
mexico_mun <- project(vect(file.path(dir$raw, "Shapefiles", "inegi_municipios", "Municipios_2010_5.shp")), crs(chirps))
chirps <- crop(chirps, mexico_mun)
chirps <- mask(chirps, mexico_mun)
cat("Mask done, elapsed:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n") 

# Check rainfall raster via plot, it looks good
# Note there was a proj_identify warning but I think its okay, you checked it out a bit
nlyr(chirps)
if (interactive()) plot(chirps[[1]]) 

# Extract rainfall values inside each municipality   
# Inside each municipality, and exact = TRUE adds a "fraction" column which is the share of each cell inside the municipality.
# This is necessary as many municipalities are smaller than one 0.05 degree cell and would be dropped otherwise. Hence we need a share
rainfall_extract <- extract(chirps, mexico_mun, cells=TRUE, exact=TRUE)
cat("Extract done, elapsed:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")   


# Check
head(rainfall_extract)
dim(rainfall_extract)

# Add the 5 digit INEGI municipality key cvegeo (2 digit state + 3 digit municipality), the same key built in the MxFLS cleaning scripts  
# First we are giving indices to rows, then we are creating the 5 digit INEGI (paste glues together, sprint pads 0's in front, and integer changes them to integer)
mun_info <- data.frame(ID = 1:nrow(mexico_mun), cvegeo = paste0(sprintf("%02d", as.integer(mexico_mun$CVE_ENT)), sprintf("%03d", as.integer(mexico_mun$CVE_MUN))))
rainfall_extract <- merge(rainfall_extract, mun_info, by="ID")
head(rainfall_extract)


# Extract longitude and latitude, only need one "stack" because the location doesn't change. Also the first column refers to long and the second to latitude.
coords <- xyFromCell(chirps[[1]], rainfall_extract$cell)
rainfall_extract$lon <- coords[,1]
rainfall_extract$lat <- coords[,2]
head(rainfall_extract[,c("cvegeo","lon","lat")])   

#===============================================================================
# 2). Convert rainfall into panel format
#===============================================================================

# Convert to data.table so data.table::melt is used and variable.name/value.name are respected
setDT(rainfall_extract)
# Reshape the columns to only be the listed ones below and put all the values in these columns. Cvegeo and fraction are used instead of the previous State
rainfall_long <- melt(rainfall_extract, id.vars=c("ID","cell","cvegeo","fraction","lon","lat"), variable.name="date", value.name="rainfall", variable.factor=FALSE)
# Convert "YYYY-MM" text to a proper Date, first of the month
rainfall_long[, date := as.Date(paste0(date, "-01"))]
cat("Reshape done, elapsed:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")   

# Check
head(rainfall_long)
dim(rainfall_long)   # 41949864 8
length(unique(rainfall_long$cvegeo))   # 2456
range(rainfall_long$date)   # Should be 1981-01-01 to 2013-12-01

#===============================================================================
# 3). Cleaning rainfall to be municipality-year
#===============================================================================

# Split up the date between year and month in order to average
rainfall_long[, `:=`(year = year(date), month = month(date))]

# Average across all grid cells within a municipality, weighting each cell by the share of it inside the municipality. Tells us how much it rained in a municipality
mun_month <- as.data.frame(rainfall_long[, .(rainfall = weighted.mean(rainfall, fraction, na.rm = TRUE)), by = .(cvegeo, year, month)]) 
cat("Municipality-month done, elapsed:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n") 

# Sum the 12 monthly values to get total annual rainfall per municipality-year   
# Full year is used because Mexico's natural disasters could fall across the whole year,
mun_year <- aggregate(rainfall ~ cvegeo + year, data = mun_month, FUN = sum, na.rm = TRUE)   
names(mun_year)[names(mun_year) == "rainfall"] <- "annual_rainfall"   


# Baseline 1981-2004 mean and SD of annual rainfall, per municipality, meaning each municipality compared to its own history
# Ends in 2004 because the DID shock window starts at the MxFLS-2 interviews 2005, so the baseline never overlaps the years being tested. 
# The 20th percentile is added for drought according to Manisha Shah and Bryce Millett Steinberg
# CAN CHECK THIS RAINFALL DROUGHT SHOCK DEFINITION
mun_climatology <- aggregate(annual_rainfall ~ cvegeo, data = mun_year[mun_year$year <= 2004, ], FUN = function(x) c(mean = mean(x), sd = sd(x), p20 = unname(quantile(x, 0.20))))  
# We fix and rename the columns
mun_climatology <- do.call(data.frame, mun_climatology)  
names(mun_climatology) <- c("cvegeo", "long_run_mean", "long_run_sd", "long_run_p20")   

# Merge climatology back onto the municipality-year panel. This is a many to one merge   
mun_year <- merge(mun_year, mun_climatology, by = "cvegeo")   
# Drop the 1981-2000 years, so the exported panel covers only 2001-2013.
mun_year <- mun_year[mun_year$year >= 2001, ] 
mun_year <- mun_year[order(mun_year$cvegeo, mun_year$year), ] 

# Similar to Ponnusamy. We will build the three rainfall variables
mun_year$log_rainfall <- log(mun_year$annual_rainfall)   
mun_year$positive_shock <- as.integer(mun_year$annual_rainfall >= mun_year$long_run_mean + mun_year$long_run_sd)   
mun_year$negative_shock <- as.integer(mun_year$annual_rainfall <= mun_year$long_run_mean - mun_year$long_run_sd)   
# Standardized deviation from the municipality's own baseline, "normalized deviation", allows +/-2 SD extreme categories in Stata 
mun_year$rain_z <- (mun_year$annual_rainfall - mun_year$long_run_mean) / mun_year$long_run_sd   

# Add binary drought shock for the DID. Drought = annual rainfall below the 20th percentile of the municipality's own 1981-2004 rainfall by Shah & Steinberg
# 20th percentile of each district's own history. 
# We code treat2009 = drought in calendar year 2009 in the child's 2005 municipality (cvegeo05). Note that the 2005-2008 droughts in that same municipality 
# enter as controls (drought2005-drought2008 x post).
mun_year$drought <- as.integer(mun_year$annual_rainfall < mun_year$long_run_p20)   
attr(mun_year$drought, "label") <- "drought (indicator)"
row.names(mun_year) <- NULL   
# Add municipality name next to cvegeo so we know which municipality it is. so match() looks up the INEGI name, NOM_MUN, for each cvegeo, iconv changed in order to fix accents
mun_year$mun_name <- iconv(mexico_mun$NOM_MUN[match(mun_year$cvegeo, mun_info$cvegeo)], from = "latin1", to = "UTF-8")  
mun_year <- mun_year[, c("cvegeo", "mun_name", setdiff(names(mun_year), c("cvegeo", "mun_name")))]   

# Export into stata
write_dta(mun_year, file.path(dir$clean, "Rainfall", "municipality_year_rainfall.dta")) 
cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")   


