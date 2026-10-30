*******************************************************************
* Subanalyses: DID by rural vs urban 
* Project: MS Thesis
*******************************************************************

* Same DID, triple difference, event studies and robustness checks as DID_Regs.do, run separately for rural and urban children
* Rural is the child's 2005 household locality "rural_05", below 2,500 people, the same rural variable DID_Regs.do uses as a 2002 control

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
global clean_mxfls   "$root/Data/Cleaned/MxFLS"
global clean_rain    "$root/Data/Cleaned/Rainfall"
global out_regs      "$root/Output/Regression_Estimates"
global out_figs      "$root/Output/Figures"
global out_logs      "$root/Output/Logs"
foreach d in "$root/Output" "$out_regs" "$out_figs" "$out_logs" {
    capture mkdir "`d'"
}
foreach d in "$out_regs" "$out_logs" "$out_figs" {
    foreach s in Main Appendix {
        capture mkdir "`d'/`s'"
    }
}
capture mkdir "$out_figs/Miscellaneous"

* If the master already opened the log, the capture fails quietly and the master keeps logging. If this file is run on its own, it opens the same log here and closes it at the end, so the DID and the event study are always in the log.
capture log using "$out_logs/Appendix/subanalysis_rural_did_regs.log", replace text name(sub_rural_did)
local own_log = (_rc == 0)

*******************************************************************
* Load DID analysis data
*******************************************************************

* Loading in data
use "$clean_mxfls/mxfls_did_analysis.dta", clear

* Rural from the child's 2005 household, built the same way as the 2005 controls in Build_DID_Analysis_Data.do 
gen tmp05 = rural if year == 2005
bysort pid_link_canonical: egen rural_05 = max(tmp05)
drop tmp05
label variable rural_05 "Child's 2005 household is rural, below 2,500 people, 1=Rural 0=Urban"

* Event study leads and lags: 2002 = -2 is the lead, 2005 = -1 is the reference, 2009 = 0
gen lead_2_treat = treat2009 * (year == 2002)
gen lag_0_treat  = treat2009 * (year == 2009)
gen lead_2_remit = remittance_hh_02 * (year == 2002)
gen lag_0_remit  = remittance_hh_02 * (year == 2009)
gen lead_2_ddd   = treat2009 * remittance_hh_02 * (year == 2002)
gen lag_0_ddd    = treat2009 * remittance_hh_02 * (year == 2009)

* IMPORTANT NOTE IS THIS IS JUST A ROBUSTNESS CHECK, IT DROPS THE SAMPLE SIZE SIGNIFICANTLY 

* We run the same event study as above but only children measured in all three waves with the outcome nonmissing in each. This way the 2002 lead and the 2009 coefficient come from the same children
foreach v in haz waz bmiz {
    bysort pid_link_canonical: egen n3_`v' = total(e_`v' == 1 & !missing(`v'))
    gen e_`v'_all3 = (e_`v' == 1 & !missing(`v') & n3_`v' == 3)
}

*******************************************************************
* Fixed Effects & Controls
*******************************************************************

* Same as DID_Regs do file, except rural_02_post is taken out of controls02 because rural does not vary within a subsample
global FE_month "pid_link_canonical int_ym"
global did_controls "drought2005_post drought2006_post drought2007_post drought2008_post"
global controls02 "dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post"
global controls05 "child_male_05_post mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post"
global controls05_birth "father_height_cm_05_post father_edu_years_05_post father_age_at_birth_05_post father_employed_05_post birth_order_05_post oldest_child_05_post birth_weight_kg_05_post delivered_home_05_post predominant_breastfed_05_post iron_pregnancy_05_post premature_05_post"

* reghdfe and its dependencies, same as DID_Regs.do
capture which reghdfe
if _rc ssc install reghdfe
capture which ftools
if _rc ssc install ftools
capture which require
if _rc ssc install require

*******************************************************************
* Regressions by rural
*******************************************************************

* Everything below runs once for rural and once for urban children. 
* Note that: Figures and event study files get _rural or _urban at the end of the name
foreach s in rural urban {
if "`s'" == "rural" local rural "rural_05 == 1"
if "`s'" == "urban" local rural "rural_05 == 0"
di as result _n "==================== `s' ===================="

*******************************************************************
* DID regressions: 2009 drought
*******************************************************************

* Now we do the regressions again splitting on rural 
foreach v in haz waz bmiz {
    * Model 1: Basic DID
    reghdfe `v' treat_post if s_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    * Model 2: DID controlling for the other droughts in the window
    reghdfe `v' treat_post $did_controls if s_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    * Model 3: DID with baseline 2002 household controls and 2005 child and mother controls
    reghdfe `v' treat_post $did_controls $controls02 $controls05 if s_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
	* Model 4: Model 3 with the father and birth controls, smaller sample
    reghdfe `v' treat_post $did_controls $controls02 $controls05 $controls05_birth if s_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    * Model 5: Triple difference, does 2002 remittance receipt buffer the drought
    reghdfe `v' treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    * Total drought effect for children in 2002 remittance households
    lincom treat_post + treat_post_remit02
}

*******************************************************************
* Parallel trends and event study
*******************************************************************

* Parallel trends graph and event study split on rural
foreach v in haz waz bmiz {
    preserve
    keep if e_`v' == 1 & `rural'
    collapse (mean) `v', by(treat2009 year)
    twoway ///
        (connected `v' year if treat2009 == 0, mc(navy) lc(navy) ms(O)) ///
        (connected `v' year if treat2009 == 1, mc(cranberry) lc(cranberry) ms(T)), ///
        xline(2009, lp(dash) lc(black)) ///
        xlabel(2002 2005 2009) ///
        legend(order(1 "Control" 2 "Treated (2009 drought)") pos(6) row(1)) ///
        xtitle("Wave") ytitle("Mean `v'") ///
        title("Mean `v' by wave, treated vs control, `s'")
    graph export "$out_figs/Miscellaneous/trends_`v'_`s'.png", replace
    restore
}

* Event-study TWFE DID
foreach v in haz waz bmiz {
    reghdfe `v' lead_2_treat lag_0_treat if e_`v' == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    test lead_2_treat

    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Appendix/es_`v'_`s'.dta", replace
    post handle (-2) (_b[lead_2_treat]) (_se[lead_2_treat]) (_b[lead_2_treat] - 1.96*_se[lead_2_treat]) (_b[lead_2_treat] + 1.96*_se[lead_2_treat])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_treat]) (_se[lag_0_treat]) (_b[lag_0_treat] - 1.96*_se[lag_0_treat]) (_b[lag_0_treat] + 1.96*_se[lag_0_treat])
    postclose handle

    use "$out_regs/Appendix/es_`v'_`s'.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("Event study: `v', `s'") legend(off)
    graph export "$out_figs/Miscellaneous/event_study_`v'_`s'.png", replace
    restore
}

* c. Triple-difference event study
foreach v in haz waz bmiz {
    reghdfe `v' lead_2_ddd lag_0_ddd lead_2_treat lag_0_treat lead_2_remit lag_0_remit if e_`v' == 1 & !missing(remittance_hh_02) & `rural', absorb($FE_month) vce(cluster cvegeo05)
    test lead_2_ddd
    lincom lag_0_treat + lag_0_ddd

    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Appendix/es_ddd_`v'_`s'.dta", replace
    post handle (-2) (_b[lead_2_ddd]) (_se[lead_2_ddd]) (_b[lead_2_ddd] - 1.96*_se[lead_2_ddd]) (_b[lead_2_ddd] + 1.96*_se[lead_2_ddd])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_ddd]) (_se[lag_0_ddd]) (_b[lag_0_ddd] - 1.96*_se[lag_0_ddd]) (_b[lag_0_ddd] + 1.96*_se[lag_0_ddd])
    postclose handle

    use "$out_regs/Appendix/es_ddd_`v'_`s'.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("DDD event study (drought x 2002 remittances): `v', `s'") legend(off)
    graph export "$out_figs/Miscellaneous/event_study_ddd_`v'_`s'.png", replace
    restore
}

*******************************************************************
* Robustness tests
*******************************************************************

*******************************************************************
* a. Children measured in all three waves
*******************************************************************

foreach v in haz waz bmiz {
    reghdfe `v' treat_post if s_`v' == 1 & e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    reghdfe `v' treat_post $did_controls if s_`v' == 1 & e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    reghdfe `v' treat_post $did_controls $controls02 $controls05 if s_`v' == 1 & e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    reghdfe `v' treat_post $did_controls $controls02 $controls05 $controls05_birth if s_`v' == 1 & e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    reghdfe `v' treat_post remittance_hh_02_post treat_post_remit02 $did_controls $controls02 $controls05 if s_`v' == 1 & e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    lincom treat_post + treat_post_remit02
}

foreach v in haz waz bmiz {
    preserve
    keep if e_`v'_all3 == 1 & `rural'
    collapse (mean) `v', by(treat2009 year)
    twoway ///
        (connected `v' year if treat2009 == 0, mc(navy) lc(navy) ms(O)) ///
        (connected `v' year if treat2009 == 1, mc(cranberry) lc(cranberry) ms(T)), ///
        xline(2009, lp(dash) lc(black)) ///
        xlabel(2002 2005 2009) ///
        legend(order(1 "Control" 2 "Treated (2009 drought)") pos(6) row(1)) ///
        xtitle("Wave") ytitle("Mean `v'") ///
        title("Mean `v' by wave, children in all 3 waves, `s'")
    graph export "$out_figs/Miscellaneous/trends_`v'_all3_`s'.png", replace
    restore

    reghdfe `v' lead_2_treat lag_0_treat if e_`v'_all3 == 1 & `rural', absorb($FE_month) vce(cluster cvegeo05)
    test lead_2_treat

    preserve
    postfile handle int etime double b se lb ub using "$out_regs/Appendix/es_`v'_all3_`s'.dta", replace
    post handle (-2) (_b[lead_2_treat]) (_se[lead_2_treat]) (_b[lead_2_treat] - 1.96*_se[lead_2_treat]) (_b[lead_2_treat] + 1.96*_se[lead_2_treat])
    post handle (-1) (0) (0) (0) (0)
    post handle (0) (_b[lag_0_treat]) (_se[lag_0_treat]) (_b[lag_0_treat] - 1.96*_se[lag_0_treat]) (_b[lag_0_treat] + 1.96*_se[lag_0_treat])
    postclose handle

    use "$out_regs/Appendix/es_`v'_all3_`s'.dta", clear
    sort etime
    twoway ///
        (rcap ub lb etime, lcolor(black)) ///
        (scatter b etime, mcolor(blue) msymbol(O)), ///
        xlabel(-2 "2002" -1 "2005 (ref)" 0 "2009") ///
        xline(-1, lpattern(dash) lcolor(red)) ///
        yline(0, lpattern(dash) lcolor(black)) ///
        xtitle("Event time (wave)") ytitle("Coefficient") ///
        title("Event study, children in all 3 waves: `v', `s'") legend(off)
    graph export "$out_figs/Miscellaneous/event_study_`v'_all3_`s'.png", replace
    restore
}
}

if `own_log' log close sub_rural_did
