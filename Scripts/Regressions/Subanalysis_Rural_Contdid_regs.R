#===============================================================================
# MS Thesis
# Subanalyses: Continuous DID by rural vs urban 
#===============================================================================

# Same continuous DID as Contdid_regs.R, run separately for rural and urban children
# Rural is the child's 2005 household locality, below 2,500 people
# Sections follow Contdid_regs.R. Case 3 in section 5 is not rerun because we aren't using it
# Functions come from functions.R

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 1). Preliminary -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Clean up workspace
rm(list=ls())
while (sink.number() > 0) sink()
start_time <- Sys.time()

# Load packages
options(repos = c(CRAN = "https://cloud.r-project.org"))
want <- c("haven","contdid","did","sandwich","BMisc")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need, dependencies = TRUE)
if (packageVersion("did") < "2.5.0") install.packages("did")
if (packageVersion("contdid") < "0.1.1") {
  if (packageVersion("ptetools") < "1.0.1") install.packages("ptetools", type = "source")
  if (!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")
  remotes::install_github("bcallaway11/contdid", upgrade = "never")
}
if (!"fwildclusterboot" %in% installed.packages()[,"Package"]) {
  if (!"remotes" %in% installed.packages()[,"Package"]) install.packages("remotes")
  remotes::install_version("summclust", version = "0.7.2", upgrade = "never")
  remotes::install_version("fwildclusterboot", version = "0.13.0", upgrade = "never")
}
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
dir.create(dir$output_regs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir$output_figure, recursive = TRUE, showWarnings = FALSE)
dir$output_logs <- paste(dir$root,"/Output/Logs",sep="")
for (d in c(dir$output_regs, dir$output_figure, dir$output_logs)) for (s in c("Main","Appendix")) dir.create(paste(d,"/",s,sep=""), recursive = TRUE, showWarnings = FALSE)
for (d in c(dir$output_figure, dir$output_logs)) dir.create(paste(d,"/Miscellaneous",sep=""), recursive = TRUE, showWarnings = FALSE)

# Save to this log
sink(paste(dir$output_logs,"/Appendix/subanalysis_rural_contdid_regs.log",sep=""), split = TRUE)

# Same settings as Contdid_regs.R but we also add rural
misc <- list()
misc$outcomes <- c("haz","waz","bmiz")
misc$spline_knots <- 1
misc$spline_degree <- 3
misc$biters <- 1000
misc$age_cut <- c(haz = 16, waz = 11, bmiz = 16)
misc$rural <- c(rural = 1, urban = 0)

# Load the functions from functions.R
source(paste(dir$root,"/Scripts/Regressions/functions.R",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Load the panel, the same data as Contdid_regs.R  -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# First we are just reading in the data, need less variables because don't need controls, add rural here
panel <- read_dta(paste(dir$clean,"/MxFLS/mxfls_did_analysis.dta",sep=""),
                  col_select = c("pid_link_canonical","year","cvegeo05","age","haz","waz","bmiz","remittance_hh_02",
                                 "treat2009","treat_early","rain_z2009","rain_deficit2009","any_deficit2009","rural"))
panel <- as.data.frame(panel)

# Same sample rules as DID_Regs.do. A 2005 municipality with 2009 rainfall (treat2009 not missing), and we drop every child measured before November 2009 (treat_early = 1)
panel <- panel[!is.na(panel$treat2009) & panel$treat_early == 0, ]

# Child's rural status from the 2005 row, the same as rural_05 in Stata. Coded in original as 1 = rural, 0 = urban
# We keep the 2005 row, one row per child, then attach it to every row of the child 
m05 <- panel[panel$year == 2005 & as.numeric(panel$rural) %in% c(0, 1), c("pid_link_canonical","rural")]
m05$rural_05 <- as.integer(as.numeric(m05$rural) == 1)
panel <- merge(panel, m05[, c("pid_link_canonical","rural_05")], by = "pid_link_canonical")

# Dose and any deficit
panel$D <- panel$rain_deficit2009
panel$D_any <- as.integer(panel$any_deficit2009)
rain09 <- unique(panel[, c("cvegeo05","rain_z2009")])

# ---- Build the panels for rural and urban ----

# This code builds the same two-wave and three-wave panels as Contdid_regs.R for each outcome, once for rural and once for urban children
panel_all <- panel
df_list <- list()
df3_list <- list()
for (s in names(misc$rural)) {
  panel <- panel_all[panel_all$rural_05 == misc$rural[[s]], ]
  for (v in misc$outcomes) {
    df_list[[paste0(v, "_", s)]] <- make_df(v, c(2005, 2009))
    df3_list[[paste0(v, "_", s)]] <- make_df(v, c(2002, 2005, 2009))
  }
}
panel <- panel_all

# Check
sapply(df_list, count_children)
sapply(df3_list, count_children)
elapsed("Data")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 3). Case 1: Dose aggregation, ATT(d) and ACRT(d), by rural -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Assumption here: the dose is D = max(0, -z), so any 2009 rain below the 1981-2004 mean z < 0 counts as treated, this includes even small deficits that are normal years which is unrealistic

cd_res <- list()
for (s in names(misc$rural)) {
  for (v in misc$outcomes) {
    k <- paste0(v, "_", s)
    df <- df_list[[k]]
    set.seed(1234)
    cd_res[[k]] <- cont_did(
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
    cat("\n==================== Case 1, dose aggregation:", v, "|", s, "====================\n")
    print(summary(cd_res[[k]]))
    save_plot(ggcont_did(cd_res[[k]], type = "att"), paste0("Miscellaneous/contdid_att_", k, ".png"))
    save_plot(ggcont_did(cd_res[[k]], type = "acrt"), paste0("Miscellaneous/contdid_acrt_", k, ".png"))
  }
}
elapsed("Case 1")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 4). Case 2: Event study aggregations, by rural -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This is for the event study. Thus we need children in all three waves.
# Lists to store the result outcomes. For parallel trends and strict parallel trends
cd_res_es_level <- list()
cd_res_es_slope <- list()
for (s in names(misc$rural)) {
  for (v in misc$outcomes) {
    k <- paste0(v, "_", s)
    df3 <- df3_list[[k]]

    # Event study aggregation for ATT, tests parallel trends
    set.seed(1234)
    cd_res_es_level[[k]] <- cont_did(
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
    cd_res_es_slope[[k]] <- cont_did(
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

    cat("\n==================== Case 2, ATT event study:", v, "|", s, "====================\n")
    print(summary(cd_res_es_level[[k]]))
    cat("\n==================== Case 2, ACRT event study:", v, "|", s, "====================\n")
    print(summary(cd_res_es_slope[[k]]))
    save_plot(ggcont_did(cd_res_es_level[[k]]), paste0("Miscellaneous/contdid_es_att_", k, ".png"))
    save_plot(ggcont_did(cd_res_es_slope[[k]]), paste0("Miscellaneous/contdid_es_acrt_", k, ".png"))
  }
}
elapsed("Case 2")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 6). Binned scatter and municipality-clustered SEs, by rural -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# The purpose of this section is correct the package that gets the estimates right but the standard errors wrong for my setting.
# The reasoning for the incorrect errors is that we are at the child level but our SEs need to be at our treatment level which is at municipality level. AKA contdid gets the estimates right but its standard errors treat every child as independent.
# We use Callaways code on his github
# Additionally this creates the Callaway's Binscatter

results <- data.frame()

cd_res_cluster <- list()
for (s in names(misc$rural)) {
  for (v in misc$outcomes) {
    k <- paste0(v, "_", s)
    df <- df_list[[k]]
    fd <- first_diff(df)
    trt <- fd[fd$D > 0, ]
    unt <- fd[fd$D == 0, ]

    cat("\n==================== Clustered estimates:", v, "|", s, "====================\n")
    cat("Untreated (D = 0):", nrow(unt), "children in", length(unique(unt$cvegeo05)), "municipalities |",
        "Treated (D > 0):", nrow(trt), "children in", length(unique(trt$cvegeo05)), "municipalities\n")
    print(round(quantile(trt$D, c(0, .1, .25, .5, .75, .9, 1)), 2))

    # ---- a. Binned scatter ----
    trt$bin <- cut(trt$D, breaks = unique(quantile(trt$D, seq(0, 1, .1))), include.lowest = TRUE)
    bins <- aggregate(cbind(D, dY) ~ bin, data = trt, FUN = mean)
    bins$n <- as.numeric(table(trt$bin)[as.character(bins$bin)])
    png(paste(dir$output_figure,"/Miscellaneous/contdid_binscatter_",k,".png",sep=""), width = 6.5, height = 4.5, units = "in", res = 200)
    plot(bins$D, bins$dY, pch = 21, bg = "grey70", cex = 1 + 2 * bins$n / max(bins$n),
         xlab = "2009 rainfall deficit (SDs below the 1981-2004 mean)", ylab = paste0("Change in ", v, ", 2005 to 2009"),
         main = paste0("Change in ", v, " by 2009 rainfall deficit, ", s), ylim = range(c(bins$dY, mean(unt$dY))))
    abline(h = mean(unt$dY), lty = 2)
    legend("bottomleft", legend = c("Treated, deciles of the deficit", "Untreated mean (D = 0)"),
           pch = c(21, NA), pt.bg = c("grey70", NA), lty = c(NA, 2), bty = "n", cex = 0.8)
    dev.off()

    # ---- b. and c. Clustered SEs ----
    res <- cd_res[[k]]
    n <- nrow(fd)
    n1 <- nrow(trt)
    b_att.d <- mboot_mun(res$att.d_inffunc, df)
    b_acrt.d <- mboot_mun(res$acrt.d_inffunc, df)
    b_att <- mboot_mun(res$overall_att_inffunc, df)
    b_acrt <- mboot_mun(res$overall_acrt_inffunc * n / n1, df)
    cd_res_cluster[[k]] <- list(att.d_se = b_att.d$se, att.d_crit.val = b_att.d$crit.val,
                                acrt.d_se = b_acrt.d$se, acrt.d_crit.val = b_acrt.d$crit.val,
                                overall_att_se = b_att$se, overall_acrt_se = b_acrt$se)

    cat("Average ATT  =", round(res$overall_att, 3), "(", round(b_att$se, 3), ")   [contdid SE:", round(res$overall_att_se, 3), "]\n")
    cat("Average ACRT =", round(res$overall_acrt, 3), "(", round(b_acrt$se, 3), ")   [contdid SE:", round(res$overall_acrt_se, 3), "]\n")
    add_result(v, paste0("2005-2009, ", s), "Average ATT (parallel trends)", res$overall_att, res$overall_att_se, b_att$se, nrow(fd))
    add_result(v, paste0("2005-2009, ", s), "Average ACRT (needs strong parallel trends)", res$overall_acrt, res$overall_acrt_se, b_acrt$se, nrow(fd))

    # ---- d. Plots ----
    xlab <- "2009 rainfall deficit (SDs below the 1981-2004 mean)"
    plot_bands(res$dose, res$att.d, b_att.d$se, b_att.d$crit.val, paste0("Miscellaneous/contdid_att_clustered_", k, ".png"),
               xlab, "ATT", paste0("ATT(d), ", v, ", ", s, ": 95% bands clustered by municipality"))
    plot_bands(res$dose, res$acrt.d, b_acrt.d$se, b_acrt.d$crit.val, paste0("Miscellaneous/contdid_acrt_clustered_", k, ".png"),
               xlab, "ACRT", paste0("ACRT(d), ", v, ", ", s, ": 95% bands clustered by municipality"))
  }
}

elapsed("Clustered section")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 7). Remittance buffer: Case 1 separately by 2002 remittance receipt, by rural -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Since contdid has no covariates, we run Case 1 separately for children in 2002 remittance and non-remittance households, with parallel trends assumed within each group
# The clustered difference in the average ATT between the two groups is the continuous-DID version of the Stata triple difference
# Essentially part one of the loop calculates this drought effect separately for each group. And if you subtract it you would get the second part of the loop

# Tests whether 2002 remittances protected children from the 2009 drought. It loops over the three outcomes HAZ, WAZ, BMIZ, and for each one does two things.
# 1: Runs CGBS continuous DID, splits children into non-remittance and remittance households, then estimates the drought effect aka the average ATT within each group with a municipality-clustered SE.
# 2: Puts both groups back together and estimates the difference between those two effects, remittance minus non-remittance, also clustered. If it’s positive and significant for HAZ, the drought hurt remittance children less.


cd_res_remit <- list()
for (s in names(misc$rural)) {
  for (v in misc$outcomes) {
    k <- paste0(v, "_", s)

    # ---- 1. Each group separately ----
    for (r in c(0, 1)) {
      name <- paste0(k, "_remit", r)
      df_r <- df_list[[k]]
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

      fd_r <- first_diff(df_r)
      fit_r <- lm(dY ~ D_any, data = fd_r)
      b_r <- coef(fit_r)[2]
      se_r <- sqrt(sandwich::vcovCL(fit_r, cluster = ~cvegeo05, type = "HC1")[2, 2])
      cat("\n==================== Remittance buffer:", v, "|", s, "| 2002 remittances =", r, "====================\n")
      cat("Untreated:", sum(fd_r$D == 0), "children in", length(unique(fd_r$cvegeo05[fd_r$D == 0])), "municipalities |",
          "Treated:", sum(fd_r$D > 0), "children\n")
      print(summary(cd_res_remit[[name]]))
      cat("Clustered average ATT =", round(b_r, 3), "(", round(se_r, 3), ")\n")
      add_result(v, paste0("2005-2009, ", s), paste0("Average ATT, 2002 remittances = ", r), b_r, cd_res_remit[[name]]$overall_att_se, se_r, nrow(fd_r))
    }

    # ---- 2. Difference between the groups ----
    fd_b <- first_diff(df_list[[k]][!is.na(df_list[[k]]$remittance_hh_02), ])
    fit_int <- lm(dY ~ D_any * remittance_hh_02, data = fd_b)
    b_int <- coef(fit_int)["D_any:remittance_hh_02"]
    V_int <- sandwich::vcovCL(fit_int, cluster = ~cvegeo05, type = "HC1")
    se_int <- sqrt(V_int["D_any:remittance_hh_02", "D_any:remittance_hh_02"])
    cat("Difference in average ATT (remittance minus non-remittance) =", round(b_int, 3), "(", round(se_int, 3), ")\n")
    add_result(v, paste0("2005-2009, ", s), "Average ATT difference: remittance minus non-remittance", b_int, NA, se_int, nrow(fd_b))
  }
}
elapsed("Remittance section")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 8). SPI drought classes, with wild cluster bootstrap p-values, by rural -----
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

spi_results <- data.frame()

cd_res_spi <- list()
for (s in names(misc$rural)) {
  for (v in misc$outcomes) {
    k <- paste0(v, "_", s)
    df_spi <- add_spi(df_list[[k]])
    fd_spi <- first_diff(df_spi)
    near_normal <- function(fd) fd$D_spi == 0
    cat("\n==================== SPI classes:", v, "|", s, "====================\n")
    cat("Near normal (untreated):", sum(fd_spi$D_spi == 0), "children in", length(unique(fd_spi$cvegeo05[fd_spi$D_spi == 0])), "municipalities |",
        "Dry (z <= -1):", sum(fd_spi$dry), "children in", length(unique(fd_spi$cvegeo05[fd_spi$dry == 1])), "municipalities\n")

    # ---- a. contdid with the SPI dose ----
    d_pos <- df_spi$D_spi[df_spi$D_spi > 0]
    set.seed(1234)
    cd_res_spi[[k]] <- cont_did(
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
    print(summary(cd_res_spi[[k]]))
    save_plot(ggcont_did(cd_res_spi[[k]], type = "att"), paste0("Miscellaneous/contdid_spi_att_", k, ".png"))

    res_spi <- cd_res_spi[[k]]
    n_trt <- sum(fd_spi$D_spi > 0)
    b_att_spi.d <- mboot_mun(res_spi$att.d_inffunc, df_spi)
    b_att_spi <- mboot_mun(res_spi$overall_att_inffunc, df_spi)
    b_acrt_spi <- mboot_mun(res_spi$overall_acrt_inffunc * nrow(fd_spi) / n_trt, df_spi)
    cat("a. Clustered: Average ATT =", round(res_spi$overall_att, 3), "(", round(b_att_spi$se, 3), ") |",
        "Average ACRT =", round(res_spi$overall_acrt, 3), "(", round(b_acrt_spi$se, 3), ")\n")
    row <- data.frame(outcome = v, specification = paste0("contdid with the SPI dose (a.), ", s), term = "Average ACRT (needs strong parallel trends)",
                      coef = round(res_spi$overall_acrt, 4), se_clustered = round(b_acrt_spi$se, 4),
                      p_clustered = round(2 * pnorm(-abs(res_spi$overall_acrt / b_acrt_spi$se)), 3), p_wild_bootstrap = NA,
                      n_children = nrow(fd_spi), children_in_group = n_trt,
                      municipalities_in_group = length(unique(fd_spi$cvegeo05[fd_spi$D_spi > 0])),
                      untreated_municipalities = length(unique(fd_spi$cvegeo05[fd_spi$D_spi == 0])), row.names = NULL)
    spi_results <- rbind(spi_results, row)

    plot_bands(res_spi$dose, res_spi$att.d, b_att_spi.d$se, b_att_spi.d$crit.val, paste0("Miscellaneous/contdid_spi_att_clustered_", k, ".png"),
               "2009 rainfall deficit (SDs below the 1981-2004 mean), z <= -1", "ATT",
               paste0("ATT(d), SPI dose, ", v, ", ", s, ": 95% bands clustered by municipality"))

    # ---- b. Dry vs near normal ----
    cat("b. Average ATT, dry vs near normal:\n")
    spi_reg(fd_spi, "dry", v, paste0("Average ATT: dry (z <= -1) vs near normal, ", s), near_normal)

    # ---- c. SPI groups ----
    cat("c. SPI groups vs near normal:\n")
    spi_reg(fd_spi, c("mod_dry","sev_dry"), v, paste0("Groups: moderately dry, severely dry or worse, ", s), near_normal)
    cat("Extremely dry shown alone (few municipalities, SE not reliable):\n")
    spi_reg(fd_spi, c("mod_dry","sev_only","ext_dry"), v, paste0("Groups: all three SPI dry classes, ", s), near_normal)

    # ---- f. Remittance buffer ----
    cat("f. Remittance buffer, dry x 2002 remittances:\n")
    fd_rem <- fd_spi[!is.na(fd_spi$remittance_hh_02), ]
    spi_reg(fd_rem, c("dry","remittance_hh_02","dry:remittance_hh_02"), v, paste0("Remittance buffer: dry x 2002 remittances, ", s), near_normal)

    # ---- g. Main average ATT ----
    cat("g. Wild bootstrap check of the main average ATT (any deficit vs none):\n")
    no_deficit <- function(fd) fd$D == 0
    spi_reg(first_diff(df_list[[k]]), "D_any", v, paste0("Main average ATT from 6): any deficit vs none, ", s), no_deficit)
  }
}
elapsed("SPI section")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 9). Save the headline results -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

cat("\n==================== Summary ====================\n")
print(results, row.names = FALSE)
write.csv(results, paste(dir$output_regs,"/Appendix/subanalysis_rural_contdid_results.csv",sep=""), row.names = FALSE)
cat("\n==================== SPI classes and wild bootstrap ====================\n")
print(spi_results, row.names = FALSE)
write.csv(spi_results, paste(dir$output_regs,"/Appendix/subanalysis_rural_contdid_spi_results.csv",sep=""), row.names = FALSE)

cat("Total runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")
sink()

# The end
