#===============================================================================
# MS Thesis
# Continuous DID functions
#===============================================================================

# Contdid_regs.R loads this file in section 1. This file should not be run on its own

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Build the panel the same way as DID_Regs.do -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Create a balanced panel for one outcome and set of waves, in order to use the format: id, G, D, time_period, Y
make_df <- function(v, waves) {
  # Keep rows that are in waves 2005 and 2009, age is under the limit, and is not missing
  df <- panel[panel$year %in% waves & panel$age < misc$age_cut[[v]] & !is.na(panel[[v]]), ]
  # This code below is a safeguard for dups but we don't need it
  df <- df[!duplicated(df[, c("pid_link_canonical","year")]), ]
  # Counts how many rows a child has with the table function
  n <- table(df$pid_link_canonical)
  # Only keeping children in both waves
  df <- df[df$pid_link_canonical %in% names(n)[n == length(waves)], ]
  # We give the child an id aka an integer value
  df$id <- as.integer(factor(df$pid_link_canonical))
  # Turns years e.g. 2002, 2005, 2009, into period numbers 1, 2, 3
  df$time_period <- match(df$year, sort(waves))
  # Treated children are first treated in the last period, never treated children get G = 0
  df$G <- ifelse(df$D > 0, length(waves), 0)
  # Copies outcomes into column called Y
  df$Y <- df[[v]]
  # This code reorders by child and time period and then we keep the appropriate variables for this df
  keep <- c("id","G","D","time_period","Y","D_any","cvegeo05","remittance_hh_02")
  df <- df[order(df$id, df$time_period), keep]
  row.names(df) <- NULL
  df
}

# Dose values for ATT(d) and ACRT(d), from the smallest to the largest positive dose
# The package's default, 10th-99th percentiles of the dose, gives a distorted ATT(d) curve in contdid 0.1.1, so we use the full range of the dose instead
# Why: cont_did_acrt() fits the spline with splines2::bSpline(dose[dose > 0], ...) but evaluates it with splines2::bSpline(dvals, ...) without Boundary.knots (cont_did.R lines 404 and 411),
# so splines2 puts the spline's end points at range(dvals) instead of the range of the dose. With dvals running from the smallest to the largest treated dose the two ranges are the same
dose_values <- function(df) {
  # Keeps only treated children with a positive deficit
  d_pos <- df$D[df$time_period == 1 & df$D > 0]
  # Makes 100 evenly spaced numbers from the smallest to the largest deficit
  seq(min(d_pos), max(d_pos), length.out = 100)
}

# Save contdid plot to a png
save_plot <- function(p, name) {
  png(paste(dir$output_figure,"/",name,sep=""), width = 6.5, height = 4.5, units = "in", res = 200)
  print(p)
  dev.off()
}

# Elapsed time, printed at the end of each section
elapsed <- function(section) {
  cat(section, "done, elapsed:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")
}

# Count children, untreated and treated in each panel. Used as check after the panels are built
count_children <- function(df) {
  c(children = length(unique(df$id)),
    untreated = sum(df$G == 0 & df$time_period == 1),
    treated = sum(df$G > 0 & df$time_period == 1))
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 6). Paper's binned scatter and municipality-clustered SEs -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# first_diff() turns each child’s before-and-after rows into a single row holding how much the child’s outcome changed between the two periods for example HAZ in 2009 minus HAZ in 2005
# This follows Callaway's code exactly
first_diff <- function(df) {
  df <- BMisc::make_balanced_panel(df, "id", "time_period")
  df$.dy <- BMisc::get_first_difference(df, "id", "Y", "time_period")
  maxT <- max(df[["time_period"]])
  post_data <- subset(df, df[["time_period"]] == maxT)
  # Renamed to dY, the name the rest of this script uses
  post_data$dY <- post_data$.dy
  post_data
}

# We use mboot() function because the contdid doesn't support it and our SEs are wrong because they treat every child as independent whereas mboot allows us to cluster at municipality
# We revise upon his code to let us cluster at the treatment level. But we also utilize his ready made package to do this.
# Note the following options are specified in this package
mboot_mun <- function(inf_func, df) {
  # We grab the municipality of the child aka the last row, and we get these three vars
  units <- df[df$time_period == max(df$time_period), c("id", "time_period", "cvegeo05")]
  # Order it by units for alignment
  units <- units[order(units$id), ]
  # Stops if rows are misaligned
  stopifnot(nrow(units) == NROW(inf_func))
  set.seed(1234)
  # mboot function from callaway
  params <- list(data = units, idname = "id", tname = "time_period", clustervars = "cvegeo05",
                 biters = misc$biters, alp = 0.05, panel = TRUE, true_repeated_cross_sections = FALSE,
                 allow_unbalanced_panel = FALSE, faster_mode = FALSE)
  did::mboot(inf_func, DIDparams = params)
}

# Adds a row to results, the data frame created in section 6 of Contdid_regs.R
add_result <- function(outcome, sample, estimate, coef, se_contdid, se_cluster, n) {
  row <- data.frame(outcome = outcome, sample = sample, estimate = estimate, coef = round(coef, 4),
                    se_contdid_bootstrap = round(se_contdid, 4),
                    se_clustered_municipality = round(se_cluster, 4), n_children = n)
  results <<- rbind(results, row)
}

# Plot a dose curve where dark = pointwise 95% CI, light = uniform 95% band, both clustered by municipality
plot_bands <- function(dose, est, se, cv, name, xlab, ylab, main) {
  png(paste(dir$output_figure,"/",name,sep=""), width = 6.5, height = 4.5, units = "in", res = 200)
  plot(dose, est, type = "n", ylim = range(c(est - cv * se, est + cv * se)), xlab = xlab, ylab = ylab, main = main)
  polygon(c(dose, rev(dose)), c(est - cv * se, rev(est + cv * se)), col = "grey90", border = NA)     # uniform band
  polygon(c(dose, rev(dose)), c(est - 1.96 * se, rev(est + 1.96 * se)), col = "grey70", border = NA) # pointwise CI
  lines(dose, est, lwd = 2)
  abline(h = 0, lty = 2)
  legend("bottomleft", legend = c("Pointwise 95% CI", "Uniform 95% band"), fill = c("grey70", "grey90"),
         border = NA, bty = "n", cex = 0.8)
  dev.off()
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 8). SPI drought classes, with wild cluster bootstrap p-values -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Here we add the SPI columns
add_spi <- function(df) {
  # We match on the 2005 municipalities to 2009 rainfall and also just keeping not wet z-scores
  df$z <- rain09$rain_z2009[match(df$cvegeo05, rain09$cvegeo05)]
  df <- df[df$z < 1, ]
  # This code puts a z-score as positive if the drought is below -1 and 0 if it is not and then the next one tells us when treated
  df$D_spi <- ifelse(df$z <= -1, -df$z, 0)
  df$G <- ifelse(df$D_spi > 0, max(df$time_period), 0)
  # Just coding as binary
  df$dry <- as.integer(df$z <= -1)
  df$mod_dry <- as.integer(df$z <= -1 & df$z > -1.5)
  df$sev_dry <- as.integer(df$z <= -1.5)
  df$sev_only <- as.integer(df$z <= -1.5 & df$z > -2)
  df$ext_dry <- as.integer(df$z <= -2)
  df
}

# Run regression, save each coefficient with clustered and wild bootstrap p-values
spi_reg <- function(fd, terms, v, label, untreated) {
  # Reformulate() takes the list of names and turns it into the formula, very efficient, this is running a regression, then we cluster
  fit <- lm(reformulate(terms, "dY"), data = fd)
  V <- sandwich::vcovCL(fit, cluster = ~cvegeo05, type = "HC1")
  # Number of municipalities, used for degrees of freedom
  G <- length(unique(fd$cvegeo05))
  # Calculating coefficient, se, p-value, and the p-value with the bootstrap
  for (term in terms) {
    b <- coef(fit)[term]
    se <- sqrt(V[term, term])
    p_cl <- 2 * pt(-abs(b / se), G - 1)
    dqrng::dqset.seed(1234)
    bt <- fwildclusterboot::boottest(fit, param = term, clustid = "cvegeo05", B = 9999, type = "webb", impose_null = TRUE, conf_int = FALSE)
    p_wb <- fwildclusterboot::pval(bt)
    # Count how many children are in the group a coefficient belongs to, for the children_in_group and municipalities_in_group columns of the output. They don’t affect any estimate
    parts <- strsplit(term, ":")[[1]]
    in_group <- rowSums(fd[, parts, drop = FALSE] == 1) == length(parts)
    n_mun <- length(unique(fd$cvegeo05[in_group]))
    # Adding a row in the results column
    row <- data.frame(outcome = v, specification = label, term = term, coef = round(b, 4), se_clustered = round(se, 4),
                      p_clustered = round(p_cl, 3), p_wild_bootstrap = round(p_wb, 3), n_children = nrow(fd),
                      children_in_group = sum(in_group), municipalities_in_group = n_mun,
                      untreated_municipalities = length(unique(fd$cvegeo05[untreated(fd)])), row.names = NULL)
    spi_results <<- rbind(spi_results, row)
    cat(sprintf(" %-28s %7.3f (%5.3f)  p clustered = %.3f | p wild bootstrap = %.3f | %d children in %d municipalities\n",
                term, b, se, p_cl, p_wb, sum(in_group), n_mun))
  }
}


