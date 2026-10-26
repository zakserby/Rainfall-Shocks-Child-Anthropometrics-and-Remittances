*******************************************************************
* DID: Drought between MxFLS-2 and MxFLS-3 and child health
* Project: MS Thesis
*******************************************************************

clear all
set more off

*******************************************************************
* Directories
*******************************************************************
* The root folder is set once in 00_Master.do. If this file is run on its own, the root is found by going up from Stata's working folder until the folder that has Scripts/00_Master.do
if "$root" == "" {
    local r = subinstr("`c(pwd)'", "\", "/", .)
    while "`r'" != "" & !fileexists("`r'/Scripts/00_Master.do") {
        local r = substr("`r'", 1, strrpos("`r'", "/") - 1)
    }
    global root "`r'"
}
global raw_mxfls     "$root/Data/Raw/MxFLS"
global raw_who       "$root/Data/Raw/WHO_reference"
global clean_mxfls   "$root/Data/Cleaned/MxFLS"
global clean_who     "$root/Data/Cleaned/WHO_zscores"
global clean_rain    "$root/Data/Cleaned/Rainfall"
global clean_banxico "$root/Data/Cleaned/Banxico"
global out_sumstats  "$root/Output/Tables"
global out_regs      "$root/Output/Regression_Estimates"
global out_logs      "$root/Output/Logs"
global out_figs      "$root/Output/Figures"
capture mkdir "$out_figs"
* Create the folders this project writes to. Note: capture skips quietly if a folder already exists)
foreach d in "$root/Data" "$root/Data/Cleaned" "$clean_mxfls" "$clean_who" "$clean_rain" "$root/Output" "$out_sumstats" "$out_regs" "$out_logs" {
    capture mkdir "`d'"
}
foreach d in "$out_regs" "$out_logs" "$out_figs" {
    foreach s in Main Appendix {
        capture mkdir "`d'/`s'"
    }
}

* If the master already opened the log, the capture fails quietly and the master keeps logging. If this file is run on its own, it opens the same log here and closes it at the end, so the DID and the event study are always in the log.
capture log using "$out_logs/Main/did_regs.log", replace text name(regs_did)
local own_log = (_rc == 0)

*******************************************************************
* Load DID analysis data
*******************************************************************

* Loading in data
use "$clean_mxfls/mxfls_did_analysis.dta", clear


*******************************************************************
* Fixed Effects & Controls
*******************************************************************

global FE_did "pid_link_canonical year"
* Child FE and interview year-month FE, so treated and control children measured in the same month are compared. 
global FE_month "pid_link_canonical int_ym"
* Controls for other droughts. Must be interacted with post
global did_controls "drought2005_post drought2006_post drought2007_post drought2008_post"
* Baseline 2002 household controls. Must be interacted with post
global controls02 "dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post rural_02_post"
* Child and mother 2005 controls. Must be interacted with post. These are filled for about 95% of the DID children
global controls05 "child_male_05_post mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post"
* Father and Book IV birth controls 2005. Must be interacted with post. Only filled for about 65% to 80% of the DID children so we keep them separate
global controls05_birth "father_height_cm_05_post father_edu_years_05_post father_age_at_birth_05_post father_employed_05_post birth_order_05_post oldest_child_05_post birth_weight_kg_05_post delivered_home_05_post predominant_breastfed_05_post iron_pregnancy_05_post premature_05_post"

*******************************************************************
* DID regressions: 2009 drought
*******************************************************************

* reghdfe needed for absorbing fixed effects; install once if not already present
capture which reghdfe
if _rc ssc install reghdfe
capture which ftools
if _rc ssc install ftools
capture which require
if _rc ssc install require

* NOTE: Standard errors are clustered at the 2005 municipality, where treatment is assigned. STILL IN PROCESS OF THINKING ABOUT CVEOGEO INSTEAD OF CVEOGEO05. Reasoning for 2005 is because cveogeo changes for movers and thus the clusters would not stay fixed over time

*******************************************************************
* CHILD HEIGHT (haz), age < 16
*******************************************************************

* Model 1: Basic DID
reghdfe haz treat_post if s_haz == 1, absorb($FE_month) vce(cluster cvegeo05) 
* Model 2: DID controlling for the other droughts in the window
reghdfe haz treat_post $did_controls if s_haz == 1, absorb($FE_month) vce(cluster cvegeo05)
* Model 3: DID with baseline 2002 household controls and 2005 child and mother controls
reghdfe haz treat_post $did_controls $controls02 $controls05 if s_haz == 1, absorb($FE_month) vce(cluster cvegeo05) 
* Model 4: Model 4 with the father and birth controls, smaller sample
* Note: We don't use this as the main model because father data is often missing when the father has migrated - often to the US, so requiring it drops migrant households, the ones most relevant to the remittance buffer
reghdfe haz treat_post $did_controls $controls02 $controls05 $controls05_birth if s_haz == 1, absorb($FE_month) vce(cluster cvegeo05)
* Model 5: Triple difference, does 2002 remittance receipt buffer the drought
reghdfe haz treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_haz == 1, absorb($FE_month) vce(cluster cvegeo05) 


* Total drought effect for children in 2002 remittance households. the triple difference coefficient alone only tests the gap between the two groups
lincom treat_post + treat_post_remit02



*******************************************************************
* CHILD WEIGHT (waz), age < 11
*******************************************************************

reghdfe waz treat_post if s_waz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe waz treat_post $did_controls if s_waz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe waz treat_post $did_controls $controls02 $controls05 if s_waz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe waz treat_post $did_controls $controls02 $controls05 $controls05_birth if s_waz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe waz treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_waz == 1, absorb($FE_month) vce(cluster cvegeo05) 
lincom treat_post + treat_post_remit02


*******************************************************************
* CHILD BMI (bmiz), age < 16
*******************************************************************

reghdfe bmiz treat_post if s_bmiz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe bmiz treat_post $did_controls if s_bmiz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe bmiz treat_post $did_controls $controls02 $controls05 if s_bmiz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe bmiz treat_post $did_controls $controls02 $controls05 $controls05_birth if s_bmiz == 1, absorb($FE_month) vce(cluster cvegeo05) 
reghdfe bmiz treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_bmiz == 1, absorb($FE_month) vce(cluster cvegeo05) 
lincom treat_post + treat_post_remit02



*******************************************************************
* Parallel trends and event study
*******************************************************************

* Event study sample. Here we use the same restrictions as the main DID, but all three waves
* e_haz, e_waz and e_bmiz are now built in Build_DID_Analysis_Data.do

* a. Average outcome by wave for treated and control, just a simple parallel trends
foreach v in haz waz bmiz {
    preserve
    keep if e_`v' == 1
    collapse (mean) `v', by(treat2009 year)
    twoway ///
        (connected `v' year if treat2009 == 0, mc(navy) lc(navy) ms(O)) ///
        (connected `v' year if treat2009 == 1, mc(cranberry) lc(cranberry) ms(T)), ///
        xline(2009, lp(dash) lc(black)) ///
        xlabel(2002 2005 2009) ///
        legend(order(1 "Control" 2 "Treated (2009 drought)") pos(6) row(1)) ///
        xtitle("Wave") ytitle("Mean `v'") ///
        title("Mean `v' by wave, treated vs control")
    * Saved to the Figures folder
    graph export "$out_figs/Appendix/trends_`v'.png", replace
    restore
}

* b. Event-study TWFE DID. Note that Event time: 2002 = -2 is lead, 2005 = -1 is reference, 2009 = 0 is treatment period
gen lead_2_treat = treat2009 * (year == 2002)
gen lag_0_treat  = treat2009 * (year == 2009)

foreach v in haz waz bmiz {
    reghdfe `v' lead_2_treat lag_0_treat if e_`v' == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Here we conduct a pretrend test where we expect lead to be around 0
    test lead_2_treat

    * Save coefficients and CIs for plotting, the reference period is posted as 0
    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Main/es_`v'.dta", replace
    post handle (-2) (_b[lead_2_treat]) (_se[lead_2_treat]) (_b[lead_2_treat] - 1.96*_se[lead_2_treat]) (_b[lead_2_treat] + 1.96*_se[lead_2_treat])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_treat]) (_se[lag_0_treat]) (_b[lag_0_treat] - 1.96*_se[lag_0_treat]) (_b[lag_0_treat] + 1.96*_se[lag_0_treat])
    postclose handle

    * Plot the event study coefficients
    use "$out_regs/Main/es_`v'.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("Event study: `v'") legend(off)
    * Saved to the Figures folder
    graph export "$out_figs/Main/event_study_`v'.png", replace
    restore
}


* c. Triple-difference event study
* This tests whether children in 2002 remittance households in drought areas were already on a different health path before the 2009 drought. Want to see that the 2002 lead is roughly 0. If 0, then any difference after 2009 can more credibly be read as remittances buffering the drought, not pre-existing differences between remittance and non-remittance families

gen lead_2_remit = remittance_hh_02 * (year == 2002)
gen lag_0_remit  = remittance_hh_02 * (year == 2009)
gen lead_2_ddd   = treat2009 * remittance_hh_02 * (year == 2002)
gen lag_0_ddd    = treat2009 * remittance_hh_02 * (year == 2009)
label variable lead_2_ddd "2009 drought x 2002 remittances x 2002 wave (pre-trend lead)"
label variable lag_0_ddd  "2009 drought x 2002 remittances x 2009 wave (buffer effect)"

foreach v in haz waz bmiz {
    reghdfe `v' lead_2_ddd lag_0_ddd lead_2_treat lag_0_treat lead_2_remit lag_0_remit if e_`v' == 1 & !missing(remittance_hh_02), absorb($FE_month) vce(cluster cvegeo05)
    * Pre-trend test for the DDD. The lead should be about 0
    test lead_2_ddd
    lincom lag_0_treat + lag_0_ddd
    * Save the DDD coefficients and CIs for plotting
    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Appendix/es_ddd_`v'.dta", replace
    post handle (-2) (_b[lead_2_ddd]) (_se[lead_2_ddd]) (_b[lead_2_ddd] - 1.96*_se[lead_2_ddd]) (_b[lead_2_ddd] + 1.96*_se[lead_2_ddd])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_ddd]) (_se[lag_0_ddd]) (_b[lag_0_ddd] - 1.96*_se[lag_0_ddd]) (_b[lag_0_ddd] + 1.96*_se[lag_0_ddd])
    postclose handle
	* Plot the DDD event study coefficients
    use "$out_regs/Appendix/es_ddd_`v'.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("DDD event study (drought x 2002 remittances): `v'") legend(off)
    graph export "$out_figs/Appendix/event_study_ddd_`v'.png", replace
    restore
}


*******************************************************************
* Robustness tests
*******************************************************************

*******************************************************************
* a. Children measured in all three waves
*******************************************************************

* IMPORTANT NOTE IS THIS IS JUST A ROBUSTNESS CHECK, IT DROPS THE SAMPLE SIZE SIGNIFICANTLY 

* We run the same event study as above but only children measured in all three waves with the outcome nonmissing in each. This way the 2002 lead and the 2009 coefficient come from the same children
foreach v in haz waz bmiz {
	* Here we count each childs rows to make sure they are in all three waves. This is done by generating a flag column called "all3" for children seen in all 3 waves. First we count the total that are not missing in each row, then we flag the row if not missing and the child present in all 3 rows.
    bysort pid_link_canonical: egen n3_`v' = total(e_`v' == 1 & !missing(`v'))
    gen e_`v'_all3 = (e_`v' == 1 & !missing(`v') & n3_`v' == 3)
    label variable e_`v'_all3 "Event study sample, children in all 3 waves, `v'"
}


* Now we do the regressions again looking at just kids in all three waves
foreach v in haz waz bmiz {
    * Model 1: Basic DID
    reghdfe `v' treat_post if s_`v' == 1 & e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Model 2: DID controlling for the other droughts in the window
    reghdfe `v' treat_post $did_controls if s_`v' == 1 & e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Model 3: DID with baseline 2002 household controls and 2005 child and mother controls
    reghdfe `v' treat_post $did_controls $controls02 $controls05 if s_`v' == 1 & e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Model 4: Model 3 with the father and birth controls, smaller sample
    reghdfe `v' treat_post $did_controls $controls02 $controls05 $controls05_birth if s_`v' == 1 & e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Model 5: Triple difference, does 2002 remittance receipt buffer the drought
    reghdfe `v' treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_`v' == 1 & e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05)
    * Total drought effect for children in 2002 remittance households
    lincom treat_post + treat_post_remit02
}

* Parallel trends graph and event study for the children in all 3 waves
foreach v in haz waz bmiz {
    preserve
    keep if e_`v'_all3 == 1
    collapse (mean) `v', by(treat2009 year)
    twoway ///
        (connected `v' year if treat2009 == 0, mc(navy) lc(navy) ms(O)) ///
        (connected `v' year if treat2009 == 1, mc(cranberry) lc(cranberry) ms(T)), ///
        xline(2009, lp(dash) lc(black)) ///
        xlabel(2002 2005 2009) ///
        legend(order(1 "Control" 2 "Treated (2009 drought)") pos(6) row(1)) ///
        xtitle("Wave") ytitle("Mean `v'") ///
        title("Mean `v' by wave, children in all 3 waves")
    graph export "$out_figs/Appendix/trends_`v'_all3.png", replace
    restore

    reghdfe `v' lead_2_treat lag_0_treat if e_`v'_all3 == 1, absorb($FE_month) vce(cluster cvegeo05) 
    * Pre-trend test: lead should be about 0
    test lead_2_treat

    * Save coefficients and CIs for plotting, the reference period is posted as 0
    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Appendix/es_`v'_all3.dta", replace
    post handle (-2) (_b[lead_2_treat]) (_se[lead_2_treat]) (_b[lead_2_treat] - 1.96*_se[lead_2_treat]) (_b[lead_2_treat] + 1.96*_se[lead_2_treat])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_treat]) (_se[lag_0_treat]) (_b[lag_0_treat] - 1.96*_se[lag_0_treat]) (_b[lag_0_treat] + 1.96*_se[lag_0_treat])
    postclose handle

    * Plot the event-study coefficients
    use "$out_regs/Appendix/es_`v'_all3.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("Event study, children in all 3 waves: `v'") legend(off)
    * Saved to the Figures folder
    graph export "$out_figs/Appendix/event_study_`v'_all3.png", replace
    restore
}



*******************************************************************
* b. Propensity score matching, currently on hold
*******************************************************************

* NOTE: We are considering combining propensity score matching with the DID, so that children in 2002 remittance and non-remittance households look the same on observable pre-drought characteristics 
* NOTE: We tried it using logit propensity score, psmatch2 nearest neighbor with replacement and a caliper of 0.05, matching separately in drought and control municipalities, then the triple difference weighted by _weight, but the result was not that nice.



if `own_log' log close regs_did






