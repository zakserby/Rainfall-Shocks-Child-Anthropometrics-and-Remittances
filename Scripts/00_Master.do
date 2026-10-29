*******************************************************************
* MASTER FILE: runs the project's scripts in order
* Project: MS Thesis (Mexico: natural disasters, remittances, child health)
*******************************************************************
clear all
set more off
capture log close _all

*******************************************************************
* 1). Directories
*******************************************************************
* Every script builds its paths from $root. It is found automatically by going up from Stata's working folder until the folder that has Scripts/00_Master.do, so nothing needs to change on a new computer
* Note: Double-clicking 00_Master.do sets Stata's working folder to Scripts. If Stata was already open, cd to the project or Scripts folder first
local r = subinstr("`c(pwd)'", "\", "/", .)
while "`r'" != "" & !fileexists("`r'/Scripts/00_Master.do") {
    local r = substr("`r'", 1, strrpos("`r'", "/") - 1)
}
global root "`r'"
if "$root" == "" {
    di as error "Project folder not found. cd to the project or Scripts folder and rerun, or set global root to the project folder here"
    exit 601
}

global scripts       "$root/Scripts"
global raw_mxfls     "$root/Data/Raw/MxFLS"
global raw_who       "$root/Data/Raw/WHO_reference"
global clean_mxfls   "$root/Data/Cleaned/MxFLS"
global clean_who     "$root/Data/Cleaned/WHO_zscores"
global clean_rain    "$root/Data/Cleaned/Rainfall"
global out_sumstats  "$root/Output/Tables"
global out_regs      "$root/Output/Regression_Estimates"
global out_logs      "$root/Output/Logs"

* Create the output folders if they don't exist yet. Note: capture skips quietly if they do
foreach d in "$root/Data" "$root/Data/Cleaned" "$clean_mxfls" "$clean_who" "$clean_rain" "$root/Output" "$out_sumstats" "$out_regs" "$out_logs" "$root/Output/Figures" {
    capture mkdir "`d'"
}
* This specifies where the results, logs, stats, and figures should go 
foreach d in "$out_sumstats" "$out_regs" "$out_logs" "$root/Output/Figures" {
    foreach s in Main Appendix {
        capture mkdir "`d'/`s'"
    }
}
capture mkdir "$out_logs/Miscellaneous"
capture mkdir "$root/Output/Figures/Miscellaneous"

*******************************************************************
* 2). Select what to run, where 1 runs it and 0 doesn't
*******************************************************************

global run_R           0   // Rebuild the municipality rainfall data in R. Note: it is slow because it reads the 7.7 GB CHIRPS file). Also its output is already in Data/Cleaned 
global run_cleaning    1   // MxFLS cleans all waves and appends them 
global run_build_did   1   // Builds the DID analysis data once. Note: The summary stats, DID, continuous DID and subanalyses all load it
global run_sumstats_did 1  // DID summary statistics 
global run_map         1   // Map of the 2009 drought by municipality in R
global run_did         1   // DID regressions
global run_contdid     1   // Continuous DID in R. Slow bc the wild cluster bootstrap uses 9,999 draws
global run_sub_gender_did     1   // Subanalyses: DID regressions for boys and girls separately
global run_sub_gender_contdid 1   // Subanalyses: continuous DID for boys and girls separately. Slow for the same reason as run_contdid
global run_sub_rural_did     1   // Subanalyses: DID regressions for rural and urban children separately
global run_sub_rural_contdid 1   // Subanalyses: continuous DID for rural and urban children separately. Slow for the same reason as run_contdid
global run_regstats    1   // Exports the regression tables to LaTeX. Note: runs last because it needs the continuous DID results from R
global run_paper       1   // Puts all main and appendix tables and the main figures into Paper/Results.tex and compiles Paper/Results.pdf

* If the municipality rainfall data is not in Data/Cleaned yet, e.g. a fresh copy with only the raw data, the R rainfall step runs even if run_R is 0
capture confirm file "$clean_rain/municipality_year_rainfall.dta"
if _rc global run_R 1


* Path to Rscript (related to MacOS shell or terminal)
global Rscript "/usr/local/bin/Rscript"
* Note: This code essentially checks if Rscript is not at the path above, try the other usual Mac paths, and stop with a message if an R step is on but R is not found, so the R steps don't fail silently
capture confirm file "$Rscript"
if _rc {
    foreach p in "/opt/homebrew/bin/Rscript" "/Library/Frameworks/R.framework/Resources/bin/Rscript" {
        capture confirm file "`p'"
        if !_rc global Rscript "`p'"
    }
}
capture confirm file "$Rscript"
if _rc & ($run_R | $run_map | $run_contdid | $run_sub_gender_contdid | $run_sub_rural_contdid | $run_paper) {
    di as error "Rscript was not found. Set global Rscript to the path of Rscript on this computer"
    exit 601
}

*******************************************************************
* 3). R cleaning: municipality rainfall and drought (CHIRPS) 
*******************************************************************


* We run the rainfall data if turned on, set a directory, use shell to run the Rscript outside of R, create a log file and then 2>&1 places those errors and outputs in the log file.
if $run_R {
    cd "$scripts" 
    shell "$Rscript" "$scripts/Cleaning/Cleaned_rainfall_municipality_drought.R" > "$out_logs/Miscellaneous/Cleaned_rainfall_municipality_drought.Rout" 2>&1 
}

*********************************************************************************************************
* 4). MxFLS cleaning, Note, must be sequential as 2009 appends, then the DID analysis data is built
*********************************************************************************************************

* Similar thought process as above but now for MxFLS and for stata
if $run_cleaning {
    log using "$out_logs/Miscellaneous/cleaning_mxfls_2002.log", replace text name(clean02)
    do "$scripts/Cleaning/Cleaning Script for MxFLS_2002_03.do"
    log close clean02

    log using "$out_logs/Miscellaneous/cleaning_mxfls_2005.log", replace text name(clean05)
    do "$scripts/Cleaning/Cleaning Script for MxFLS_2005_06.do"
    log close clean05

    log using "$out_logs/Miscellaneous/cleaning_mxfls_2009.log", replace text name(clean09)
    do "$scripts/Cleaning/Cleaning Script for MxFLS_2009_12.do"
    log close clean09
}

* Builds the DID analysis data
if $run_build_did {
    log using "$out_logs/Miscellaneous/build_did_analysis_data.log", replace text name(build_did)
    do "$scripts/Cleaning/Build_DID_Analysis_Data.do"
    log close build_did
}

*******************************************************************
* 5). Summary statistics
*******************************************************************

* DID summary statistics
if $run_sumstats_did {
    log using "$out_logs/Main/summary_stats_did.log", replace text name(sumstats_did)
    do "$scripts/Statistics/Summary Stats DID.do"
    log close sumstats_did
}

* Map of the 2009 drought by municipality, it loads the DID analysis data for the sample municipalities
if $run_map {
    cd "$scripts"
    shell "$Rscript" "$scripts/Statistics/Drought_map.R" > "$out_logs/Miscellaneous/Drought_map.Rout" 2>&1
}


*******************************************************************
* 6). DID, DDD and event study regressions 
*******************************************************************

* DID regressions
if $run_did {
    log using "$out_logs/Main/did_regs.log", replace text name(regs_did)
    do "$scripts/Regressions/DID_Regs.do" 
    log close regs_did
}

*******************************************************************
* 7). Continuous DID in R 
*******************************************************************

if $run_contdid {
    cd "$scripts"
    shell "$Rscript" "$scripts/Regressions/Contdid_regs.R" > "$out_logs/Miscellaneous/Contdid_regs.Rout" 2>&1 
}

*******************************************************************
* 8). Subanalyses by child sex: boys vs girls
*******************************************************************

* Same as sections 6 and 7 but for boys and girls separately. They load the same DID analysis data, so they run last 
if $run_sub_gender_did {
    log using "$out_logs/Appendix/subanalysis_gender_did_regs.log", replace text name(sub_gender_did)
    do "$scripts/Regressions/Subanalysis_Gender_DID_Regs.do"
    log close sub_gender_did
}

* In general, as with this code below, Stata runs that R script in the background through the Mac terminal. Everything R prints, including any errors, is saved to a log file. This is if gender_contdid is switched on
if $run_sub_gender_contdid {
    cd "$scripts"
    shell "$Rscript" "$scripts/Regressions/Subanalysis_Gender_Contdid_regs.R" > "$out_logs/Miscellaneous/Subanalysis_Gender_Contdid_regs.Rout" 2>&1
}

*******************************************************************
* 9). Subanalyses by locality: rural vs urban
*******************************************************************

* Same as sections 6 and 7 but for rural and urban children separately. They load the same DID analysis data, so they run last 
if $run_sub_rural_did {
    log using "$out_logs/Appendix/subanalysis_rural_did_regs.log", replace text name(sub_rural_did) 
    do "$scripts/Regressions/Subanalysis_Rural_DID_Regs.do"
    log close sub_rural_did
}

* In general, as with this code below, Stata runs that R script in the background through the Mac terminal. Everything R prints, including any errors, is saved to a log file. This is if rural contdid is switched on
if $run_sub_rural_contdid {
    cd "$scripts"
    shell "$Rscript" "$scripts/Regressions/Subanalysis_Rural_Contdid_regs.R" > "$out_logs/Miscellaneous/Subanalysis_Rural_Contdid_regs.Rout" 2>&1 
}

*******************************************************************
* 10). Regression tables in LaTeX
*******************************************************************

* Exports the DID, continuous DID and subanalysis results to LaTeX tables in Tables/Main and Tables/Appendix. This is the same way the summary stats do file does. It reruns the Stata regressions and reads the continuous DID csv files from R, so it runs last
if $run_regstats {
    log using "$out_logs/Main/regression_stats.log", replace text name(regstats)
    do "$scripts/Statistics/Regression Stats.do"
    log close regstats
}

*******************************************************************
* 11). Results document, LaTeX and PDF
*******************************************************************

* We compile and put everything into one results document in LaTeX and PDF
if $run_paper {
    cd "$scripts"
    shell "$Rscript" "$scripts/Paper/Build_Results.R" > "$out_logs/Miscellaneous/Build_Results.Rout" 2>&1
}

di as result "Master file finished."








