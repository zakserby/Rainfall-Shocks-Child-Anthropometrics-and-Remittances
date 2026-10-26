*******************************************************************
* Build the DID analysis data: 2009 drought, controls and samples
* Project: MS Thesis
*******************************************************************

* The summary stats, DID_Regs.do, the subanalyses and the continuous DID in R all load this file it saves, mxfls_did_analysis.dta
* Note: runs after the rainfall script and the MxFLS cleaning, it needs municipality_year_rainfall.dta and mxfls_panel_2002_2005_2009.dta
* NOTE: THIS CODE CREATES THE MAIN DATASET THAT IS USED FOR BOTH DID CONTINUOUS AND BINARY DID REGRESSIONS

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

*******************************************************************
* Load MxFLS panel
*******************************************************************

use "$clean_mxfls/mxfls_panel_2002_2005_2009.dta", clear
* Rename to "year" to match the rainfall panel's merge key
rename wave year

*******************************************************************
* Merge rainfall panel, municipality-year
*******************************************************************

merge m:1 cvegeo year using "$clean_rain/municipality_year_rainfall.dta", keep(master match)
* Puts the municipality name - mun_name, from municipality_year_rainfall.dta - right after cvegeo
order mun_name, after(cvegeo)
*tab _merge
*tab year _merge
drop _merge

*******************************************************************
* Interview year
*******************************************************************

* Why interview year and not the wave label: about 18% of the 2005 wave was interviewed in 2006-07 and 38% of the 2009 wave in 2010-13, so the drought window between the waves is household specific
* Interview year = year of the last dated visit in the Book C cover file (c_conpor), the same rule the cleaning scripts use for age_days
* The 2002, 2005 and 2009 cleaning scripts all save int_date, int_year and int_month, so we don't need to rebuild them
* int_date is now the Book S measurement date
replace int_year = year if missing(int_year)


* Note: This allows us to compare treated and control kids measured in the same month. Seasonal and timing differences are taken out.
gen int_ym = ym(int_year, int_month)
format int_ym %tm
* tab year int_year

*******************************************************************
* Baseline municipality 2005 in order to avoid migration during 2009
*******************************************************************

* This code is the first part of the treatment code 

* We grab the child's 2005 municipality
gen cvegeo05 = cvegeo if year == 2005
* This code duplicates the child's 2005 municipality across rows
bysort pid_link_canonical (cvegeo05): replace cvegeo05 = cvegeo05[_N]
* We do the same procedure for municipality name in order to know the name 
gen mun_name05 = mun_name if year == 2005
bysort pid_link_canonical (mun_name05): replace mun_name05 = mun_name05[_N]
order cvegeo05 mun_name05, after(mun_name)
label variable cvegeo05 "Child's 2005 municipality"
label variable mun_name05 "Name of the child's 2005 municipality"


*******************************************************************
* Merge municipality rainfall wide on the 2005 municipality
*******************************************************************


* One row per municipality with each year's drought as its own column, merged on the child's 2005 municipality
* We use the chirps rainfall data that we cleaned, we use this range to capture the drought in 2009 and droughts beforehand in case we want to use these pre droughts as controls. 
preserve
    use "$clean_rain/municipality_year_rainfall.dta", clear
    keep cvegeo year drought rain_z
    local drought_lab : variable label drought
    keep if inrange(year, 2002, 2009)
    reshape wide drought rain_z, i(cvegeo) j(year)
    rename cvegeo cvegeo05
    tempfile rain_wide
    save `rain_wide'
restore
merge m:1 cvegeo05 using `rain_wide', keep(master match)
* tab year _merge if cvegeo05 != ""
drop _merge
foreach y in 2002 2003 2004 2005 2006 2007 2008 2009 {
    label variable drought`y' "`y' `drought_lab'"
}

*******************************************************************
* Creating treatment of 2009 drought, single shock
*******************************************************************

* We create the treatment now. 1 means the child's 2005 municipality had a drought in 2009 aka treated.
gen treat2009 = drought2009 if cvegeo05 != ""
label variable treat2009 "`: variable label drought2009'"


* D = 2009 rainfall deficit in SDs. It plays the same role as treat2009 and is time invariant. We take -rain_z2009 so a deficit is positive, then max with 0 so wet years are 0
* Note: max() in Stata ignores missing values, so we keep it missing when rain_z2009 is missing
gen double rain_deficit2009 = max(0, -rain_z2009) if !missing(rain_z2009) & cvegeo05 != ""
gen any_deficit2009 = (rain_deficit2009 > 0) if !missing(rain_deficit2009)
label variable rain_deficit2009 "2009 rainfall deficit (SD)"
label variable any_deficit2009 "Any 2009 rainfall deficit (indicator)"


* NOTE: DECIDING IF I WANT TO KEEP THIS DEFINITION OF THIS CODE CHUNK BELOW.
* NOTE: This cutoff is not addressed below in the post = 0 but rather in the gen s_haz etc.. with the filter
global cutoff = mdy(11, 1, 2009)
gen early09 = (int_date < $cutoff) if year == 2009 & !missing(int_date)
* All children including treated and control are measured before November 2009 are dropped
bysort pid_link_canonical: egen treat_early = max(early09 == 1)
label variable early09 "2009 measurement before the cutoff (November 2009)"
label variable treat_early "Child (treated or control) measured before the cutoff in 2009, dropped"

*******************************************************************
* Cleaning DID variables and controls
*******************************************************************

* Define post period and treated post
gen post = 0
replace post = 1 if int_year > 2008
gen treat_post = treat2009 * post
label variable post "MxFLS-3 (after the 2009 drought)"
label variable treat_post "2009 drought x post"

* We also control for other droughts between the 2005 and 2009 interviews, so treat_post only picks up the 2009 drought. Many children had these too. 
foreach y in 2005 2006 2007 2008 {
    gen drought`y'_post = drought`y' * post
    label variable drought`y'_post "`y' drought x post"
}

* NOTE: STILL IN PROCESS OF TINKERING WITH THIS
* 2002 household controls in wave 1 measured well before the 2009 shock
* Remittances in 2002 are kept here too. Then we interact it with post in order for it not to drop out of the regression
* NOTE: We use the child's 2005 household only to find which household they belong to, then give the child that household's 2002 values on all their rows. This is only bc some kids are born in 2003 etc..
* NOTE: Must be two step process, first we must be at the household level and then go to the kid level that is bc some kids born in 2003
foreach v in remittance_hh dwell_electricity dwell_floor_good dwell_roof_good dwell_ventilation_good rooms_per_person ln_pc_expend_month rural {
	* We create a temporary variable if the year is 2002. Then we extend to all of its rows with max
    gen tmp02 = `v' if year == 2002
    bysort folio_canonical: egen hh02 = max(tmp02)
    * Keep the household's 2002 value only on the child's 2005 row, then duplicate it to all of the child's rows
    gen tmp05 = hh02 if year == 2005
    bysort pid_link_canonical: egen `v'_02 = max(tmp05)
    drop tmp02 hh02 tmp05
    gen `v'_02_post = `v'_02 * post
    label variable `v'_02 "`: variable label `v''"
}
	
	
* Derive the child's own sex from the "gender" variable, which the cleaning scripts take from the Book C roster sex (ls04), not sa01
gen child_male = 1 if gender == 1
replace child_male = 0 if gender == 3
label variable child_male "Male (indicator)" 


* Child and parent controls at baseline before the shock still. We use 2005 and not 2002 because these are the child's own values and children born after 2002 are not in 2002.
foreach v in child_male mother_height_cm mother_edu_years mother_age_at_birth mother_employed fully_vaccinated father_height_cm father_edu_years father_age_at_birth father_employed birth_order oldest_child birth_weight_kg delivered_home predominant_breastfed iron_pregnancy premature {
    * Keep the value only on the 2005 row
    gen tmp05 = `v' if year == 2005
    * Next we duplicate 2005 value to all of the child's rows
    bysort pid_link_canonical: egen `v'_05 = max(tmp05)
    drop tmp05
    gen `v'_05_post = `v'_05 * post
}


* 2002 remittances are used as the remittance measure in the triple difference
gen treat_post_remit02 = treat_post * remittance_hh_02
label variable remittance_hh_02_post "2002 remittances x post" // this is created by the last line in the loop above
label variable treat_post_remit02 "2009 drought x post x 2002 remittances"


**********************************************************************************************
* Community bank/bank branch access in 2002 (still in process of getting it) 
* would code it similar to remittances
*********************************************************************************************


* We restrict child health outcomes to these two waves 2005 and 2009, age cut offs based on WHO and MxFLS.
gen s_haz  = inlist(year, 2005, 2009) & age < 16 & treat_early == 0 & !missing(treat2009)
gen s_waz  = inlist(year, 2005, 2009) & age < 11 & treat_early == 0 & !missing(treat2009)
gen s_bmiz = inlist(year, 2005, 2009) & age < 16 & treat_early == 0 & !missing(treat2009)

* Event study sample. Here we use the same restrictions as the main DID, but all three waves
gen e_haz  = inlist(year, 2002, 2005, 2009) & age < 16 & treat_early == 0 & !missing(treat2009)
gen e_waz  = inlist(year, 2002, 2005, 2009) & age < 11 & treat_early == 0 & !missing(treat2009)
gen e_bmiz = inlist(year, 2002, 2005, 2009) & age < 16 & treat_early == 0 & !missing(treat2009)

* NOTE: THIS IS ONLY USED BY THE SUMMARY STATS. 
* These are the children reghdfe and contdid use. The samples. LOOSELY DEFINED
foreach v in haz waz bmiz {
    gen tmp = s_`v' == 1 & !missing(`v')
    bysort pid_link_canonical: egen n_`v' = total(tmp)
    gen reg_`v' = (tmp == 1 & n_`v' == 2)
    drop tmp n_`v'
    label variable reg_`v' "Estimation sample, `v' in both 2005 and 2009"
}


*******************************************************************
* Save
*******************************************************************

compress
save "$clean_mxfls/mxfls_did_analysis.dta", replace





