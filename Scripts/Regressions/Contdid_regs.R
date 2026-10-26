#===============================================================================
# MS Thesis
# Continuous DID
#===============================================================================

# Continuous DID of the 2009 rainfall deficit on child health with the contdid package.
# Paper: Callaway, Goodman-Bacon & Sant'Anna 2024, "Difference-in-Differences with a Continuous Treatment"
# Package website: https://bcallaway11.github.io/contdid/ version 0.1.1.
# Note: Sections 3 to 5 follow the website's Case 1, 2 and 3 code.
# Sections 6 to 8 are additions for our setting: the paper's binned scatter, standard errors clustered by municipality, the remittance buffer, and the SPI drought classes
# Section 6 takes contdid's own estimates and influence functions and clusters their bootstrap by municipality with Callaway's did::mboot()

# How we define the data (id, G, D, time_period, Y)
#   id = child. Y = haz, waz or bmiz
#   D = 2009 rainfall deficit in the child's 2005 municipality, in SDs below its 1981-2004 mean. D = 0 if 2009 rain was at or above the mean
#   time_period = 1 for 2005 and 2 for 2009. Note: event study: 1 = 2002, 2 = 2005, 3 = 2009
#   G = period the child is first treated: In the two period design G for treated is 2 if D > 0. In the event study design it is 3 if D > 0, and 0 if D = 0 that means never treated
#   Same sample rules as DID_Regs.do: 2005 municipality, drop children measured before November 2009, age < 16 for HAZ, BMIZ or < 11 for WAZ

# NOTE: The package does not support covariates so no controls or interview year-month FE unlike DID_Regs.do and also unbalanced panels
# Assumptions
#   Random Sampling
#   Continuous or Multi-Valued Discrete Treatment
#   No-Anticipation and Observed Outcomes
#   Parallel trends identifies the level effects, ATT(d) and the average ATT
#   Strong parallel trends is also needed to read the slope effects, ACRT(d) and the average ACRT, or any comparison across doses as causal

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 1). Preliminary -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Clean up workspace
rm(list=ls())
# Close a log left open by an earlier run that stopped partway, technically don't need it but it's convenient to keep
while (sink.number() > 0) sink()
# Run time
start_time <- Sys.time()

# Load packages. haven reads the Stata files, fwildclusterboot runs the wild cluster bootstrap in section 8
# did, BMisc and sandwich are used in 6) to 8). did and BMisc are Callaway's packages: did::mboot() clusters the bootstrap and BMisc::get_first_difference() builds Delta Y. sandwich::vcovCL clusters the regression SEs that are not in contdid
options(repos = c(CRAN = "https://cloud.r-project.org"))
want <- c("haven","contdid","did","sandwich","BMisc")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need, dependencies = TRUE)
# Making sure we don't use an older DID
if (packageVersion("did") < "2.5.0") install.packages("did")
if (packageVersion("contdid") < "0.1.1") {
  if (packageVersion("ptetools") < "1.0.1") install.packages("ptetools", type = "source")
  if (!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")
  remotes::install_github("bcallaway11/contdid", upgrade = "never")
}

# fwildclusterboot and its dependency summclust were archived by CRAN, so they are installed once from the CRAN archive.
# On a Mac this needs the Xcode command line tools "xcode-select --install" and gfortran-12.2-universal.pkg from https://mac.r-project.org/tools/
if (!"fwildclusterboot" %in% installed.packages()[,"Package"]) {
  if (!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")
  remotes::install_version("summclust", version = "0.7.2", upgrade = "never")
  remotes::install_version("fwildclusterboot", version = "0.13.0", upgrade = "never")
}

# Load everything once and stop if any package didn't load
loaded <- sapply(c(want, "fwildclusterboot"), function(i) require(i, character.only=TRUE))
if (!all(loaded)) stop("These packages did not load: ", paste(names(loaded)[!loaded], collapse = ", "))
if (packageVersion("did") < "2.5.0") stop("did is ", packageVersion("did"), ", but 2.5.0 or newer is needed. Update it, restart R, and rerun")
rm(want, need, loaded)

# Directories
dir <- list()
dir$root <- dirname(getwd())
# If the project root is not one folder up from the working folder, go up from the working folder until the folder that has Scripts/00_Master.do
if (!dir.exists(paste(dir$root,"/Data",sep=""))) {
  dir$root <- getwd()
  while (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep="")) && dirname(dir$root) != dir$root) dir$root <- dirname(dir$root)
  if (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep=""))) stop("Project folder not found. Set the working directory to the Scripts folder (open MS Thesis.Rproj) or run 00_Master.do")
}
dir$clean <- paste(dir$root,"/Data/Cleaned",sep="")
dir$output_regs <- paste(dir$root,"/Output/Regression_Estimates",sep="")
dir$output_figure <- paste(dir$root,"/Output/Figures",sep="")
# Create the folders this script writes to (does nothing if they already exist)
dir.create(dir$output_regs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir$output_figure, recursive = TRUE, showWarnings = FALSE)
dir$output_logs <- paste(dir$root,"/Output/Logs",sep="")
for (d in c(dir$output_regs, dir$output_figure, dir$output_logs)) for (s in c("Main","Appendix")) dir.create(paste(d,"/",s,sep=""), recursive = TRUE, showWarnings = FALSE)
for (d in c(dir$output_figure, dir$output_logs)) dir.create(paste(d,"/Miscellaneous",sep=""), recursive = TRUE, showWarnings = FALSE)

# Everything printed below is also saved to this log
sink(paste(dir$output_logs,"/Main/contdid_regs.log",sep=""), split = TRUE)

# Misc just used as our list for the settings your script uses later
# We are just defining some vars here to use later. This way is more efficient
# We use the function's default of 1000 bootstrap draws for inference
misc <- list()
misc$outcomes <- c("haz","waz","bmiz")
misc$spline_knots <- 1
misc$spline_degree <- 3
misc$biters <- 1000
# Age cutoffs from DID_Regs.do
misc$age_cut <- c(haz = 16, waz = 11, bmiz = 16)

# Load the functions this script uses from functions.R, in the same folder as this script
source(paste(dir$root,"/Scripts/Regressions/functions.R",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Load the panel, the same data as DID_Regs.do  -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# First we are just reading in the data, need less variables because don't need controls
panel <- read_dta(paste(dir$clean,"/MxFLS/mxfls_did_analysis.dta",sep=""),
                  col_select = c("pid_link_canonical","year","cvegeo05","age","haz","waz","bmiz","remittance_hh_02",
                                 "treat2009","treat_early","rain_z2009","rain_deficit2009","any_deficit2009"))
panel <- as.data.frame(panel)

# Same sample rules as DID_Regs.do. A 2005 municipality with 2009 rainfall (treat2009 not missing), and we drop every child measured before November 2009 (treat_early = 1)
panel <- panel[!is.na(panel$treat2009) & panel$treat_early == 0, ]

# D = 2009 rainfall deficit in SDs. D on every row plays exactly the role treat2009 on every row plays in Stata. This is time invariant
# Note: We take - rain_z2009 in order to have a positive interpretation and then we take the max aka compare to zero and keep the greater value
panel$D <- panel$rain_deficit2009
# Any deficit, 1 if D > 0. Used for the ATT event study in 4, the remittance buffer in 7 and the wild bootstrap check in 8 g.
panel$D_any <- as.integer(panel$any_deficit2009)
rain09 <- unique(panel[, c("cvegeo05","rain_z2009")])

# ---- Functions ----

# Theses functions are in the functions r script
# make_df(), dose_values(), save_plot() and elapsed() moved to functions.R, which is loaded in section 1

# ---- Build the panels ----

# Build the two-period (2005, 2009) and three-period (2002, 2005, 2009) panels for each outcome. This is so later sections can just call upon it rather than recreating it
df_list <- list()
df3_list <- list()
for (v in misc$outcomes) {
  df_list[[v]] <- make_df(v, c(2005, 2009))
  df3_list[[v]] <- make_df(v, c(2002, 2005, 2009))
}

# Check, head(df), looks good rn
head(df_list$haz)
# count_children() is in functions.R
sapply(df_list, count_children)
sapply(df3_list, count_children)
elapsed("Data")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 3). Case 1: Dose aggregation, ATT(d) and ACRT(d) -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Assumption here: the dose is D = max(0, -z), so any 2009 rain below the 1981-2004 mean z < 0 counts as treated, this includes even small deficits that are normal years which is unrealistic

cd_res <- list()
for (v in misc$outcomes) {
  df <- df_list[[v]]
  # Data isn't random but the SE are partly random here so we set a seed, similar to callaway code
  set.seed(1234)
  cd_res[[v]] <- cont_did(
    yname = "Y",
    tname = "time_period",
    idname = "id",
    dname = "D",
    data = df,
    gname = "G",
    target_parameter = "slope",
    aggregation = "dose",
    treatment_type = "continuous",
    control_group = "notyettreated",
    biters = misc$biters,
    cband = TRUE,
    num_knots = misc$spline_knots,
    degree = misc$spline_degree,
    dvals = dose_values(df)
  )
  
  # Use this to have the results show up in the console and the log
  cat("\n==================== Case 1, dose aggregation:", v, "====================\n")
  print(summary(cd_res[[v]]))
  save_plot(ggcont_did(cd_res[[v]], type = "att"), paste0("Miscellaneous/contdid_att_", v, ".png"))
  save_plot(ggcont_did(cd_res[[v]], type = "acrt"), paste0("Miscellaneous/contdid_acrt_", v, ".png"))
}
elapsed("Case 1")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 4). Case 2: Event study aggregations -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This is for the event study. Thus we need children in all three waves.
# Lists to store the result outcomes. For parallel trends and strict parallel trends
cd_res_es_level <- list()
cd_res_es_slope <- list()
for (v in misc$outcomes) {
  df3 <- df3_list[[v]]
  
  # Event study aggregation for ATT, tests parallel trends
  set.seed(1234)
  cd_res_es_level[[v]] <- cont_did(
    yname = "Y",
    tname = "time_period",
    idname = "id",
    dname = "D_any",
    data = df3,
    gname = "G",
    target_parameter = "level",
    aggregation = "eventstudy",
    treatment_type = "continuous",
    control_group = "notyettreated",
    biters = misc$biters,
    cband = TRUE,
    num_knots = misc$spline_knots,
    degree = misc$spline_degree
  )
  
  # Event study aggregation for ACRT, tests strong parallel trends
  set.seed(1234)
  cd_res_es_slope[[v]] <- cont_did(
    yname = "Y",
    tname = "time_period",
    idname = "id",
    dname = "D",
    data = df3,
    gname = "G",
    target_parameter = "slope",
    aggregation = "eventstudy",
    treatment_type = "continuous",
    control_group = "notyettreated",
    biters = misc$biters,
    cband = TRUE,
    num_knots = misc$spline_knots,
    degree = misc$spline_degree
  )
  
  # Use this to have the results show up in the console and the log
  cat("\n==================== Case 2, ATT event study:", v, "====================\n")
  print(summary(cd_res_es_level[[v]]))
  cat("\n==================== Case 2, ACRT event study:", v, "====================\n")
  print(summary(cd_res_es_slope[[v]]))
  save_plot(ggcont_did(cd_res_es_level[[v]]), paste0("Miscellaneous/contdid_es_att_", v, ".png"))
  save_plot(ggcont_did(cd_res_es_slope[[v]]), paste0("Miscellaneous/contdid_es_acrt_", v, ".png"))
}
elapsed("Case 2")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 5). Case 3: Data-driven nonparametric model -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Note: I will likely not use this or at least not for main specification
cd_res_cck <- list()
for (v in misc$outcomes) {
  df <- df_list[[v]]
  set.seed(1234)
  cd_res_cck[[v]] <- tryCatch(cont_did(
    yname = "Y",
    tname = "time_period",
    idname = "id",
    dname = "D",
    data = df,
    gname = "G",
    target_parameter = "level",
    aggregation = "dose",
    treatment_type = "continuous",
    dose_est_method = "cck",
    control_group = "notyettreated",
    biters = misc$biters,
    cband = TRUE
  ), error = function(e) {
    cat("\nCase 3 (CCK) did not run for", v, ":", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(cd_res_cck[[v]])) {
    cat("\n==================== Case 3, CCK (not reliable here, see above):", v, "====================\n")
    print(summary(cd_res_cck[[v]]))
  }
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 6). Paper's binned scatter and municipality-clustered SEs -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# The purpose of this section is correct the package that gets the estimates right but the standard errors wrong for my setting.
# The reasoning for the incorrect errors is that we are at the child level but our SEs need to be at our treatment level which is at municipality level. AKA contdid gets the estimates right but its standard errors treat every child as independent.
# We use Callaways code on his github
# Additionally this creates the Callaway's Binscatter

# ---- Functions ----

# first_diff(), mboot_mun(), add_result() and plot_bands() moved to functions.R, which is loaded in section 1

# We are creating a data frame with the results, this inserts inputs later on. The results we want to add are defined by the function and the variables in it.
# Then we combine these results into a dataframe.
results <- data.frame()

# ---- Loop over outcomes ----

# We first are making 4 separate data frames.
# df: dataframe, 2 rows per child, 2005 and 2009
# fd: first difference, 1 row per child, with their 2005→2009 change
# trt: treated, the fd rows with a 2009 deficit (D > 0)
# unt: untreated, the fd rows with no deficit (D = 0)
cd_res_cluster <- list()
for (v in misc$outcomes) {
  df <- df_list[[v]]
  fd <- first_diff(df)
  trt <- fd[fd$D > 0, ]
  unt <- fd[fd$D == 0, ]
  
  cat("\n==================== Clustered estimates:", v, "====================\n")
  cat("Untreated (D = 0):", nrow(unt), "children in", length(unique(unt$cvegeo05)), "municipalities |",
      "Treated (D > 0):", nrow(trt), "children in", length(unique(trt$cvegeo05)), "municipalities\n")
  # Print how big the treated deficits are for children
  print(round(quantile(trt$D, c(0, .1, .25, .5, .75, .9, 1)), 2))
  
  # ---- a. Binned scatter ----
  # a. Binned scatter: mean change in the outcome by decile of the dose among treated, untreated mean as a dashed line
  # We build the bins by splitting treated children into 10 groups by dose at the 0%, 10%, …, 100% cutoffs
  trt$bin <- cut(trt$D, breaks = unique(quantile(trt$D, seq(0, 1, .1))), include.lowest = TRUE)
  # Groups treated children by their dose bin and takes the average dose and average change in the outcome within each bin. NOTE: Those averages are each circle’s x and y position on the plot.
  bins <- aggregate(cbind(D, dY) ~ bin, data = trt, FUN = mean)
  # Counts the number of children in that bin and then determines the size of the circle based on that
  bins$n <- as.numeric(table(trt$bin)[as.character(bins$bin)])
  # Simply construct the figure below
  png(paste(dir$output_figure,"/", ifelse(v %in% c("haz","waz"), "Main", "Appendix"), "/contdid_binscatter_",v,".png",sep=""), width = 6.5, height = 4.5, units = "in", res = 200)
  plot(bins$D, bins$dY, pch = 21, bg = "grey70", cex = 1 + 2 * bins$n / max(bins$n),
       xlab = "2009 rainfall deficit (SDs below the 1981-2004 mean)", ylab = paste0("Change in ", v, ", 2005 to 2009"),
       main = paste0("Change in ", v, " by 2009 rainfall deficit"), ylim = range(c(bins$dY, mean(unt$dY))))
  abline(h = mean(unt$dY), lty = 2)
  legend("bottomleft", legend = c("Treated, deciles of the deficit", "Untreated mean (D = 0)"),
         pch = c(21, NA), pt.bg = c("grey70", NA), lty = c(NA, 2), bty = "n", cex = 0.8)
  dev.off()
  
  # ---- b. and c. Clustered SEs ----
  # b. and c. We redo their SEs by clustering by municipality with did::mboot(). This is for ATT(d), ACRT(d), ATT^o and ACRT^o in contdid's own estimates from Case 1 in 3
  # Grabs outcome’s Case 1 results from cont_did() from above
  res <- cd_res[[v]]
  n <- nrow(fd)
  n1 <- nrow(trt)
  # Clustered SEs for all the dose curves, att, acrt...
  b_att.d <- mboot_mun(res$att.d_inffunc, df)
  b_acrt.d <- mboot_mun(res$acrt.d_inffunc, df)
  b_att <- mboot_mun(res$overall_att_inffunc, df)
  # We mimic a similar thing in callaway and scale by treated children. His is se_avg_acr <- sd(infl_avg_acr) / sqrt(n_treated). Note that sd is n-1 so minor difference
  b_acrt <- mboot_mun(res$overall_acrt_inffunc * n / n1, df)
  # Storing results
  cd_res_cluster[[v]] <- list(att.d_se = b_att.d$se, att.d_crit.val = b_att.d$crit.val,
                              acrt.d_se = b_acrt.d$se, acrt.d_crit.val = b_acrt.d$crit.val,
                              overall_att_se = b_att$se, overall_acrt_se = b_acrt$se)
  
  cat("Average ATT  =", round(res$overall_att, 3), "(", round(b_att$se, 3), ")   [contdid SE:", round(res$overall_att_se, 3), "]\n")
  cat("Average ACRT =", round(res$overall_acrt, 3), "(", round(b_acrt$se, 3), ")   [contdid SE:", round(res$overall_acrt_se, 3), "]\n")
  add_result(v, "2005-2009", "Average ATT (parallel trends)", res$overall_att, res$overall_att_se, b_att$se, nrow(fd))
  add_result(v, "2005-2009", "Average ACRT (needs strong parallel trends)", res$overall_acrt, res$overall_acrt_se, b_acrt$se, nrow(fd))
  
  # ---- d. Plots ----
  # d. Plot ATT(d) and ACRT(d) where dark = pointwise 95% CI, light = uniform 95% band, both clustered by municipality
  xlab <- "2009 rainfall deficit (SDs below the 1981-2004 mean)"
  
  # ATT(d): Use contdid's estimates from Case 1, with the clustered SEs and uniform critical value from b.
  plot_bands(res$dose, res$att.d, b_att.d$se, b_att.d$crit.val, paste0(ifelse(v == "haz", "Main", "Appendix"), "/contdid_att_clustered_", v, ".png"),
             xlab, "ATT", paste0("ATT(d), ", v, ": 95% bands clustered by municipality"))
  
  # ACRT(d): same plot for the slope curve
  plot_bands(res$dose, res$acrt.d, b_acrt.d$se, b_acrt.d$crit.val, paste0("Appendix/contdid_acrt_clustered_", v, ".png"),
             xlab, "ACRT", paste0("ACRT(d), ", v, ": 95% bands clustered by municipality"))
  
  # If have time can insert an event study where we cluster by municipality
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 7). Remittance buffer: Case 1 separately by 2002 remittance receipt -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Since contdid has no covariates, we run Case 1 separately for children in 2002 remittance and non-remittance households, with parallel trends assumed within each group
# The clustered difference in the average ATT between the two groups is the continuous-DID version of the Stata triple difference
# Essentially part one of the loop calculates this drought effect separately for each group. And if you subtract it you would get the second part of the loop

# Tests whether 2002 remittances protected children from the 2009 drought. It loops over the three outcomes HAZ, WAZ, BMIZ, and for each one does two things.
# 1: Runs CGBS continuous DID, splits children into non-remittance and remittance households, then estimates the drought effect aka the average ATT within each group with a municipality-clustered SE.
# 2: Puts both groups back together and estimates the difference between those two effects, remittance minus non-remittance, also clustered. If it’s positive and significant for HAZ, the drought hurt remittance children less.
cd_res_remit <- list()
for (v in misc$outcomes) {
  
  # ---- 1. Each group separately ----
  for (r in c(0, 1)) {
    name <- paste0(v, "_remit", r)
    df_r <- df_list[[v]]
    df_r <- df_r[!is.na(df_r$remittance_hh_02) & df_r$remittance_hh_02 == r, ]
    set.seed(1234)
    cd_res_remit[[name]] <- cont_did(
      yname = "Y",
      tname = "time_period",
      idname = "id",
      dname = "D",
      data = df_r,
      gname = "G",
      target_parameter = "slope",
      aggregation = "dose",
      treatment_type = "continuous",
      control_group = "notyettreated",
      biters = misc$biters,
      cband = TRUE,
      num_knots = misc$spline_knots,
      degree = misc$spline_degree,
      dvals = dose_values(df_r)
    )
    
    # Turns the two rows per child into one row with 2005→2009 change dY
    fd_r <- first_diff(df_r)
    # Regresses each child’s change on D_any (1 if any deficit, 0 if none). The D_any coefficient is (mean change of treated children) − (mean change of untreated children).
    # This is CGBS’s average ATT, the same number as contdid’s overall_att. We run it as a regression only so we can cluster its SE.
    # This just answers how much a drought hurt a remittance family or non-remittance family separately
    fit_r <- lm(dY ~ D_any, data = fd_r)
    b_r <- coef(fit_r)[2]
    se_r <- sqrt(sandwich::vcovCL(fit_r, cluster = ~cvegeo05, type = "HC1")[2, 2])
    cat("\n==================== Remittance buffer:", v, "| 2002 remittances =", r, "====================\n")
    cat("Untreated:", sum(fd_r$D == 0), "children in", length(unique(fd_r$cvegeo05[fd_r$D == 0])), "municipalities |",
        "Treated:", sum(fd_r$D > 0), "children\n")
    print(summary(cd_res_remit[[name]]))
    cat("Clustered average ATT =", round(b_r, 3), "(", round(se_r, 3), ")\n")
    add_result(v, "2005-2009", paste0("Average ATT, 2002 remittances = ", r), b_r, cd_res_remit[[name]]$overall_att_se, se_r, nrow(fd_r))
  }
  
  # ---- 2. Difference between the groups ----
  # Difference in the average ATT, remittance minus non-remittance, clustered by municipality
  # This just takes the difference and calculates if remittances serve as a buffer
  fd_b <- first_diff(df_list[[v]][!is.na(df_list[[v]]$remittance_hh_02), ])
  fit_int <- lm(dY ~ D_any * remittance_hh_02, data = fd_b)
  b_int <- coef(fit_int)["D_any:remittance_hh_02"]
  V_int <- sandwich::vcovCL(fit_int, cluster = ~cvegeo05, type = "HC1")
  se_int <- sqrt(V_int["D_any:remittance_hh_02", "D_any:remittance_hh_02"])
  cat("Difference in average ATT (remittance minus non-remittance) =", round(b_int, 3), "(", round(se_int, 3), ")\n")
  add_result(v, "2005-2009", "Average ATT difference: remittance minus non-remittance", b_int, NA, se_int, nrow(fd_b))
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 8). SPI drought classes, with wild cluster bootstrap p-values -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Here the groups follow the Standardized Precipitation Index classes from McKee, Doesken & Kleist 1993, these are the WMO standard cutoffs, chosen independently of our data
# near normal -0.99 to 0.99 = untreated (D = 0)
# moderately dry -1.0 to -1.49
# severely dry -1.5 to -1.99
# extremely dry -2 and less

# Assumption here: children with -1 < z < 0 are unaffected. This cutoff is more reasonable than z < 0 in section 3, since a deficit under 1 SD is a normal year

# Main sample drops children in wet municipalities (z >= 1, "moderately wet"), so the untreated group is only near normal
# NOTE: Extremely dry has only 6 municipalities, thus we merge it with severely dry aka "severely dry or worse", z <= -1.5
# NOTE: We use wild cluster bootstrap p-values with the fwildclusterboot package from Fischer & Roodman, the R version of Stata's boottest as we have small number of treated clusters, we rerun these simulations to see if we get similar findings
# NOTE: Each regression only needs parallel trends. We don't compare low dosage to higher dosages because comparing the groups to each other needs strong parallel trends
# NOTE: All of these SPI regressions have a bootstrap to double check with spi_reg()
# NOTE: The wild cluster bootstrap here tests the coefficient with the null imposed, which is the version designed for few clusters
# NOTE: Under parallel trends multi-valued discrete functions relative to untreated/omitted group is consistent. But compared with each other it needs stronger parallel trends assumption

# ---- Functions ----

# add_spi() and spi_reg() moved to functions.R, which is loaded in section 1
# spi_reg() adds its rows to this data frame
spi_results <- data.frame()

# ---- Loop over outcomes ----

# For each outcome, it builds the SPI data and prints how many children and municipalities are in the untreated and dry groups
cd_res_spi <- list()
for (v in misc$outcomes) {
  df_spi <- add_spi(df_list[[v]])
  fd_spi <- first_diff(df_spi)
  near_normal <- function(fd) fd$D_spi == 0
  cat("\n==================== SPI classes:", v, "====================\n")
  cat("Near normal (untreated):", sum(fd_spi$D_spi == 0), "children in", length(unique(fd_spi$cvegeo05[fd_spi$D_spi == 0])), "municipalities |",
      "Dry (z <= -1):", sum(fd_spi$dry), "children in", length(unique(fd_spi$cvegeo05[fd_spi$dry == 1])), "municipalities\n")
  
  # ---- a. contdid with the SPI dose ----
  # a. Now we implement the contdid Case 1 with the SPI dose and the more reasonable cutoff for droughts, Assumption: children with -1 < z < 0 are unaffected
  d_pos <- df_spi$D_spi[df_spi$D_spi > 0]
  set.seed(1234)
  cd_res_spi[[v]] <- cont_did(
    yname = "Y",
    tname = "time_period",
    idname = "id",
    dname = "D_spi",
    data = df_spi,
    gname = "G",
    target_parameter = "slope",
    aggregation = "dose",
    treatment_type = "continuous",
    control_group = "notyettreated",
    biters = misc$biters,
    cband = TRUE,
    num_knots = misc$spline_knots,
    degree = misc$spline_degree,
    dvals = seq(min(d_pos), max(d_pos), length.out = 100)
  )
  print(summary(cd_res_spi[[v]]))
  save_plot(ggcont_did(cd_res_spi[[v]], type = "att"), paste0("Miscellaneous/contdid_spi_att_", v, ".png"))
  
  # Here we add the clustering using mboot, exact same code as above
  res_spi <- cd_res_spi[[v]]
  n_trt <- sum(fd_spi$D_spi > 0)
  b_att_spi.d <- mboot_mun(res_spi$att.d_inffunc, df_spi)
  b_att_spi <- mboot_mun(res_spi$overall_att_inffunc, df_spi)
  b_acrt_spi <- mboot_mun(res_spi$overall_acrt_inffunc * nrow(fd_spi) / n_trt, df_spi)
  cat("a. Clustered: Average ATT =", round(res_spi$overall_att, 3), "(", round(b_att_spi$se, 3), ") |",
      "Average ACRT =", round(res_spi$overall_acrt, 3), "(", round(b_acrt_spi$se, 3), ")\n")
  row <- data.frame(outcome = v, specification = "contdid with the SPI dose (a.)", term = "Average ACRT (needs strong parallel trends)",
                    coef = round(res_spi$overall_acrt, 4), se_clustered = round(b_acrt_spi$se, 4),
                    p_clustered = round(2 * pnorm(-abs(res_spi$overall_acrt / b_acrt_spi$se)), 3), p_wild_bootstrap = NA,
                    n_children = nrow(fd_spi), children_in_group = n_trt,
                    municipalities_in_group = length(unique(fd_spi$cvegeo05[fd_spi$D_spi > 0])),
                    untreated_municipalities = length(unique(fd_spi$cvegeo05[fd_spi$D_spi == 0])), row.names = NULL)
  spi_results <- rbind(spi_results, row)
  
  # Graph note not necessary to keep if only want the results, this is the exact same code as above
  plot_bands(res_spi$dose, res_spi$att.d, b_att_spi.d$se, b_att_spi.d$crit.val, paste0(ifelse(v == "waz", "Miscellaneous", "Appendix"), "/contdid_spi_att_clustered_", v, ".png"),
             "2009 rainfall deficit (SDs below the 1981-2004 mean), z <= -1", "ATT",
             paste0("ATT(d), SPI dose, ", v, ": 95% bands clustered by municipality"))
  
  # ---- b. Dry vs near normal ----
  # b. The average effect of being in a dry municipality, z ≤ −1, compared with a near-normal one.
  cat("b. Average ATT, dry vs near normal:\n")
  # This is the reg function being run
  spi_reg(fd_spi, "dry", v, "Average ATT: dry (z <= -1) vs near normal", near_normal)
  
  # ---- c. SPI groups ----
  # c. SPI groups mod_dry and sev_dry vs near normal. We estimate each case in this regression. NOTE: That this combines sev and extreme for the first regression bc of low clusters and the second regression it separates them just to see
  cat("c. SPI groups vs near normal:\n")
  spi_reg(fd_spi, c("mod_dry","sev_dry"), v, "Groups: moderately dry, severely dry or worse", near_normal)
  cat("Extremely dry shown alone (6 municipalities, SE not reliable):\n")
  spi_reg(fd_spi, c("mod_dry","sev_only","ext_dry"), v, "Groups: all three SPI dry classes", near_normal)
  
  # ---- d. Pre-trends ----
  # d. Pre trends, same regression but on 2002 to 2005 change for each SPI group vs near normal, children in all three waves
  cat("d. Pre-trends, 2002 to 2005, each group vs near normal:\n")
  df3_spi <- add_spi(df3_list[[v]])
  pre_spi <- first_diff(df3_spi[df3_spi$time_period %in% c(1, 2), ])
  spi_reg(pre_spi, c("mod_dry","sev_dry"), v, "Pre-trend 2002-2005: moderately dry, severely dry or worse", near_normal)
  
  # ---- f. Remittance buffer ----
  # f. Remittance buffer with the SPI groups. Difference in the average ATT dry vs near normal between 2002 remittance and non-remittance children
  cat("f. Remittance buffer, dry x 2002 remittances:\n")
  # Starts from the SPI change data "fd_spi" and drops children whose 2002 remittance status is missing.
  fd_rem <- fd_spi[!is.na(fd_spi$remittance_hh_02), ]
  spi_reg(fd_rem, c("dry","remittance_hh_02","dry:remittance_hh_02"), v, "Remittance buffer: dry x 2002 remittances", near_normal)
  
  # ---- g. Main average ATT ----
  # g. Wild bootstrap p-value for the main average ATT from sections 3-6 which are any deficit vs none.
  cat("g. Wild bootstrap check of the main average ATT (any deficit vs none):\n")
  no_deficit <- function(fd) fd$D == 0
  spi_reg(first_diff(df_list[[v]]), "D_any", v, "Main average ATT from 6): any deficit vs none", no_deficit)
}
elapsed("SPI section")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 9). Save the headline results -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

cat("\n==================== Summary ====================\n")
print(results, row.names = FALSE)
write.csv(results, paste(dir$output_regs,"/Main/contdid_results.csv",sep=""), row.names = FALSE)
cat("\n==================== SPI classes and wild bootstrap ====================\n")
print(spi_results, row.names = FALSE)
write.csv(spi_results, paste(dir$output_regs,"/Main/contdid_spi_results.csv",sep=""), row.names = FALSE)

cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")
sink()

# The end
