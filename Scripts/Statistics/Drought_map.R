#===============================================================================
# MS Thesis
# Map of the 2009 drought by municipality
#===============================================================================

# Maps the 2009 rainfall in each municipality against its own 1981-2004 history, the same treatment the DID and continuous DID use
# Panel A: 2009 rainfall z-score grouped with the SPI classes of McKee, Doesken & Kleist 1993, the same classes as add_spi() in functions.R
# Panel B: 2009 drought (indicator), rainfall below the municipality's 1981-2004 20th percentile, the binary DID treatment treat2009
# Municipalities with children in the estimation sample (2005 municipality) are outlined in black
# NOTE: Needs municipality_year_rainfall.dta from Cleaned_rainfall_municipality_drought.R and mxfls_did_analysis.dta from Build_DID_Analysis_Data.do

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 1). Preliminary -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

rm(list=ls())

# Load packages
options(repos = c(CRAN = "https://cloud.r-project.org"))
packagelist <- c("sf","haven","RColorBrewer")

newpackages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(newpackages)) install.packages(newpackages)

sapply(packagelist, function(i) require(i, character.only=TRUE))
rm(newpackages, packagelist)

# Directories
dir <- list()
dir$root <- dirname(getwd())
# If the project root is not one folder up from the working folder, go up from the working folder until the folder that has Scripts/00_Master.do
if (!dir.exists(paste(dir$root,"/Data",sep=""))) {
  dir$root <- getwd()
  while (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep="")) && dirname(dir$root) != dir$root) dir$root <- dirname(dir$root)
  if (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep=""))) stop("Project folder not found. Set the working directory to the Scripts folder (open MS Thesis.Rproj) or run 00_Master.do")
}
dir$raw <- paste(dir$root,"/Data/Raw",sep="")
dir$clean <- paste(dir$root,"/Data/Cleaned",sep="")
dir$output_figure <- paste(dir$root,"/Output/Figures",sep="")
dir.create(paste(dir$output_figure,"/Appendix",sep=""), recursive = TRUE, showWarnings = FALSE)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Load the municipalities, rainfall and sample -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# INEGI municipalities, we keep its Lambert projection since it is already made for Mexico
mun <- st_read(paste(dir$raw,"/Shapefiles/inegi_municipios/Municipios_2010_5.shp",sep=""), quiet = TRUE)
# Same 5 digit cvegeo key as the rainfall script (2 digit state + 3 digit municipality)
mun$cvegeo <- paste0(sprintf("%02d", as.integer(mun$CVE_ENT)), sprintf("%03d", as.integer(mun$CVE_MUN)))

# 2009 rainfall z-score and drought for every municipality
rain <- as.data.frame(read_dta(paste(dir$clean,"/Rainfall/municipality_year_rainfall.dta",sep=""), col_select = c("cvegeo","year","rain_z","drought")))
rain <- rain[rain$year == 2009, c("cvegeo","rain_z","drought")]
mun <- merge(mun, rain, by = "cvegeo", all.x = TRUE)

# 2005 municipalities of the children in the HAZ estimation sample, the largest sample
panel <- as.data.frame(read_dta(paste(dir$clean,"/MxFLS/mxfls_did_analysis.dta",sep=""), col_select = c("cvegeo05","reg_haz")))
sample_mun <- unique(panel$cvegeo05[panel$reg_haz == 1 & panel$cvegeo05 != ""])
mun$in_sample <- mun$cvegeo %in% sample_mun
# Check, should be 133 municipalities like the summary stats
sum(mun$in_sample)

# State borders, we merge the municipalities of each state into one polygon
states <- aggregate(mun["CVE_ENT"], by = list(mun$CVE_ENT), FUN = function(x) x[1])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 3). Map -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# SPI classes, from extremely dry to wet. Near normal is -1 < z < 1
breaks <- c(-Inf, -2, -1.5, -1, 1, Inf)
labels <- c("Extremely dry, z <= -2", "Severely dry, -2 < z <= -1.5", "Moderately dry, -1.5 < z <= -1", "Near normal, -1 < z < 1", "Wet, z >= 1")
# Browns for the dry classes and teal for wet, from the BrBG palette. Near normal is grey so the droughts stand out
spi_colors <- c(brewer.pal(11, "BrBG")[c(2, 4, 5)], "grey90", brewer.pal(11, "BrBG")[9])
# Note: right = TRUE so z = -1 counts as moderately dry and z = -1.5 as severely dry, the same as add_spi()
mun$spi_class <- cut(mun$rain_z, breaks = breaks, labels = FALSE, right = TRUE)

# Function that draws one panel: municipalities filled by color, sample municipalities outlined, then states
draw_panel <- function(fill, title) {
  par(mar = c(0.5, 0.5, 2, 0.5))
  plot(st_geometry(mun), col = fill, border = NA, main = title, cex.main = 0.9)
  plot(st_geometry(mun[mun$in_sample, ]), col = NA, border = "black", lwd = 0.8, add = TRUE)
  plot(st_geometry(states), col = NA, border = "grey40", lwd = 0.4, add = TRUE)
}

png(paste(dir$output_figure,"/Appendix/drought_map_2009.png",sep=""), width = 6.5, height = 8, units = "in", res = 300)
layout(matrix(c(1,2), ncol = 1))

# Panel A: SPI classes
# NOTE: The wet municipalities (z >= 1) are almost all in central Mexico, so the any deficit continuous DID compares against mostly one region
draw_panel(ifelse(is.na(mun$spi_class), "white", spi_colors[mun$spi_class]), "A. 2009 rainfall, SPI classes")
legend("bottomleft", legend = labels, fill = spi_colors, border = NA, bty = "n", cex = 0.65)

# Panel B: binary drought, the DID treatment
# NOTE: The binary treatment and the SPI classes don't line up. About a fifth of the 2009 droughts have z above -1, near normal in the SPI classes, so the binary treatment mixes mild and severe deficits
draw_panel(ifelse(is.na(mun$drought), "white", ifelse(mun$drought == 1, "#A6611A", "grey90")), "B. 2009 drought, rainfall below the 1981-2004 20th percentile")
legend("bottomleft", legend = c("Drought", "No drought", "Municipality in the estimation sample"), fill = c("#A6611A", "grey90", NA),
       border = c(NA, NA, "black"), bty = "n", cex = 0.65)

dev.off()

# The end
