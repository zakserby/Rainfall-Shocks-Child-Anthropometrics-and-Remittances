*******************************************************************
* Cleaning Script
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
* Create the folders this project writes to. Note: capture skips quietly if a folder already exists
foreach d in "$root/Data" "$root/Data/Cleaned" "$clean_mxfls" "$clean_who" "$clean_rain" "$root/Output" "$out_sumstats" "$out_regs" "$out_logs" {
    capture mkdir "`d'"
}
* Define base directory path
global path "$raw_mxfls/hh09dta_all"

*******************************************************************
* Health Biomarkers & Anthropometrics
*******************************************************************
use "$path/hh09dta_bs/s_sa.dta", clear
* Note for now not differentiating between reported and measured I don't think
rename sa03       age
rename sa01       gender
rename sa04_22    height_cm // Revising, but currently not in use
rename sa05_21    weight_kg_gr
rename sa12_21    blood_pressure_systole
rename sa12_22    blood_pressure_diastole
rename sa12a_2    pulse
rename sa16_21    hemoglobin_g_dl
rename sa16b_21   cholesterol
rename sa16c_21   cholesterol_hdl_good
rename sa16d_21   glycosylated_hemoglobin
rename sa11_21    waist_cm
rename sa11a_21   hip_cm
* Child biomarkers
rename sa07_21    child_height_cm
rename sa08_21    child_weight_kg
* The 2 duplicate folio ls pairs are split records of the same person. We fill each row's missing values from the other row of the same person, so the row we keep has all the information. Only affects groups with 2 rows
ds folio ls, not
foreach v in `r(varlist)' {
    bysort folio ls (`v'): replace `v' = `v'[1] if missing(`v') & _N > 1
}
* Deduplicate, this removes 2 observations
sort folio ls sa21 sa26_1
by folio ls: keep if _n == _N

* Child sex must come from the Book C roster, ls04, as sex disagrees across other books but the roster is correct.
preserve
    use "$path/hh09dta_bc/c_ls.dta", clear
    keep folio ls ls04 ls01a
    * Same 2009 roster duplicate fix as below: drop the person who left (ls01a = 3) when two people share a line number
    bysort folio ls: drop if _N > 1 & ls01a == 3
    drop ls01a
    tempfile roster_sex
    save `roster_sex'
restore
merge 1:1 folio ls using `roster_sex', keep(master match) nogen
replace gender = ls04 if !missing(ls04)
drop ls04

*******************************************************************
* Calcualting HAZ, WAZ, BMIZ
*******************************************************************

* Grabbing the date
rename sa02_1 bday
rename sa02_2 bmonth
rename sa02_3 byear


* This code purely is for the missing "day" observations. It imputes a bday if we have both month and year, this allows us to keep roughly 200 observations for children get z-scores whereas before we dropped them and had no z-scores. 
gen byte bday_imputed = !inrange(bday, 1, 31) & inrange(bmonth, 1, 12) & inrange(byear, 1900, 2013)
replace bday = 15 if bday_imputed == 1

* This shows household visit/interview. Note that there is no day recorded, so we set day=15 (mid-month). We subtract the interview date from the birth date. Most of these are duplicates for a household visit, we keep the last household visit. Note that the int_date currrently is statas date not calender date so it looks odd but its okay
preserve
    use "$path/hh09dta_bc/c_conpor.dta", clear
    keep folio secuencia anio mes
    * Drop visits with no year/month first, 2009 has 1, so the "last visit" that actually has a date will be used
    drop if missing(anio, mes)
    sort folio secuencia
    by folio: keep if _n == _N
    gen int_date = mdy(mes, 15, 2000 + anio)
    keep folio int_date
    tempfile visitdate
    save `visitdate'
restore
merge m:1 folio using `visitdate', keep(master match) nogen

* Children are weighed and measured in Book S. In 2009 the Book S notes file has the exact date of the measurement session (ne02_11 day, ne02_12 month, ne02_13 year), so we use it first. If a household has no s_ne date we use the last dated Book S visit (s_conpor), day = 15 like Book C. 
* A few households have more than one s_ne session, we keep the last one, same as the last visit rule
preserve
    use "$path/hh09dta_bs/s_ne.dta", clear
    gen ne_date = mdy(ne02_12, ne02_11, 2000 + ne02_13) if inrange(ne02_13, 9, 13)
    collapse (max) ne_date, by(folio)
    tempfile nedate
    save `nedate'
    use "$path/hh09dta_bs/s_conpor.dta", clear
    keep folio secuencia anio mes
    drop if missing(anio, mes)
    sort folio secuencia
    by folio: keep if _n == _N
    gen meas_date = mdy(mes, 15, 2000 + anio)
    merge 1:1 folio using `nedate', nogen
    replace meas_date = ne_date if !missing(ne_date)
    keep folio meas_date
    tempfile measdate
    save `measdate'
restore
merge m:1 folio using `measdate', keep(master match) nogen

* Okay about 59 negative entries but this is less than .19% so we will just leave it for now. A couple might be because of the day 15 we chose above and some might be erroneous measurement




* IMPORTANT NOTE: THINKING ABOUT ADDING IN AGE AS A CONTROL???? 
gen age_days = int_date - mdy(bmonth, bday, byear)
gen age_months_exact = age_days / 30.4375
* CHANGED: Same rule as 2002 and 2005, for ALL interviews. Some Book S birth dates are wrong, or the Book C date is long after the measurement (2011-2013 interviews: 75% of ages disagree in 2013). If the computed age disagrees with the age reported in Book S (sa03) by more than about a year, set age_days to missing so the child gets no z-score
gen int int_year = year(int_date)
replace age_days = . if !missing(age, age_days) & !inrange(age_days/365.25 - age, -0.5, 1.5)

* The check above uses the Book C date because the Book S age sa03 is the roster age from the Book C visit for most of 2009 children. Children who pass get their age at the Book S measurement date, and int_date becomes that date, or the Book C date if there is none, so the interview month FE and the November 2009 cutoff use it too
replace age_days = meas_date - mdy(bmonth, bday, byear) if !missing(age_days, meas_date)
replace int_date = meas_date if !missing(meas_date)
drop meas_date
replace age_months_exact = age_days / 30.4375
drop int_year
sum age_months_exact age_days

* Just some variable recoding below
gen sex = 1 if gender == 1
replace sex = 2 if gender == 3
* Note: In 2009, 92% of children under 2 were not measured for length, even though the Book S questionnaire says to measure them lying down. So 2009 under-5 HAZ and BMIZ mostly cover ages 2-5, WAZ covers all ages. This is a data limitation, not a code problem
* 2009's standing/lying flag
gen measure = "h" if sa07a == 1
replace measure = "l" if sa07a == 3
gen str6 ageunit = "days"
* This refers to fluid swelling and is a variable in WHO but not MxFLS so everyone gets a no
gen oedema = "n"
* Every is given equal weight
gen sw = 1
tempfile z_u5 z_5to19 z_all


* WHO's official macros (not on SSC): https://www.who.int/tools/child-growth-standards/software
* and https://www.who.int/tools/growth-reference-data-for-5to19-years/application-tools
* Note: WHO reference files live in Data/Raw/WHO_reference and the z-score output goes to Data/Cleaned/WHO_zscores
global who_igrowup_dir "$raw_who/igrowup_stata"
global who_who2007_dir "$raw_who/who2007_stata"
adopath + "$who_igrowup_dir"
adopath + "$who_who2007_dir"
global who_output_dir "$clean_who"

preserve
    * Keep children under 5 and child must have at least one of weight
    keep if age_days >= 0 & age_days < 1857 & !missing(age_days, sex) & (!missing(child_weight_kg) | !missing(child_height_cm))
    gen str244 reflib = "$who_igrowup_dir"
    gen str244 datalib = "$who_output_dir"
    gen str30 datalab = "mxfls2009_u5"
    * Deleting last runs output so we don't have a silent error here
    capture erase "$who_output_dir/mxfls2009_u5_z_rc.dta"
    igrowup_restricted reflib datalib datalab sex age_days ageunit child_weight_kg child_height_cm measure oedema sw
    * This file also creates prevelance tables in the code provided from git and erases the data _zlen, _zwei, and _zbmi, so we just need to reload them 
	use "$who_output_dir/mxfls2009_u5_z_rc.dta", clear
    * The original ado computes flags but it doesn't implement the flags and remove implausible z-scores, here we implement it
    gen haz = _zlen if _flen == 0
    gen waz = _zwei if _fwei == 0
    gen bmiz = _zbmi if _fbmi == 0
	* Implement a check to count how many got flagged/excluded as out-of-range
    count if _flen == 1
    di "u5: `r(N)' HAZ obs flagged out of range (<-6 or >6)"
    count if _fwei == 1
    di "u5: `r(N)' WAZ obs flagged out of range (<-6 or >5)"
    count if _fbmi == 1
    di "u5: `r(N)' BMIZ obs flagged out of range (<-5 or >5)"
    keep folio ls haz waz bmiz
    save `z_u5'
restore


preserve
    * Keep children under 19 and child must have at least one of weight
    keep if age_months_exact >= 61 & age_months_exact <= 228 & !missing(age_days, sex) & (!missing(child_weight_kg) | !missing(child_height_cm))
    gen str244 reflib = "$who_who2007_dir"
    gen str244 datalib = "$who_output_dir"
    gen str30 datalab = "mxfls2009_5to19"
    capture erase "$who_output_dir/mxfls2009_5to19_z.dta"
    who2007 reflib datalib datalab sex age_days ageunit child_weight_kg child_height_cm oedema sw
	use "$who_output_dir/mxfls2009_5to19_z.dta", clear
    gen haz = _zhfa if _fhfa == 0
    gen waz = _zwfa if _fwfa == 0 & age_months_exact <= 120   // WHO 2007 has no WAZ past age 10
    gen bmiz = _zbfa if _fbfa == 0
    count if _fhfa == 1
    di "5to19: `r(N)' HAZ obs flagged out of range (<-6 or >6)"
    count if _fwfa == 1
    di "5to19: `r(N)' WAZ obs flagged out of range (<-6 or >5)"
    count if _fbfa == 1
    di "5to19: `r(N)' BMIZ obs flagged out of range (<-5 or >5)"
    keep folio ls haz waz bmiz
    save `z_5to19'
restore

preserve
    use `z_u5', clear
    append using `z_5to19'
    save `z_all'
restore
merge 1:1 folio ls using `z_all', keep(master match) nogen
label variable haz "Height-for-age z-score (HAZ)"
label variable waz "Weight-for-age z-score (WAZ)"
label variable bmiz "BMI-for-age z-score (BMIZ)"
label variable age "Age (years, not in regression)"


* Keep birth month and birth year for birth month and birth year fixed effects.
label variable bmonth "Birth month of child"
label variable byear "Birth year of child"
* Keep the household interview date
gen int_year = year(int_date)
gen int_month = month(int_date)
format int_date %td
label variable int_date "Measurement date (Book S session or last Book S visit, else last Book C visit)" 
label variable int_year "Interview year"
label variable int_month "Interview month"
keep folio ls age gender height_cm weight_kg_gr blood_pressure_systole blood_pressure_diastole pulse hemoglobin_g_dl cholesterol cholesterol_hdl_good glycosylated_hemoglobin child_weight_kg child_height_cm waist_cm hip_cm age_days age_months_exact haz waz bmiz bday_imputed bmonth byear int_date int_year int_month

tempfile temp_biomarkers
save `temp_biomarkers'


*******************************************************************
* Child Health Module
*******************************************************************
* NOTE LET'S HOLD OFF RUNNING THE REGRESSION WITH THIS FOR NOW
use "$path/hh09dta_b5/v_esn.dta", clear
rename esn01     health_of_child
rename esn07_1b  has_asthma
rename esn07_1d  other_chronic_illness
rename esn07_1e  obese
rename esn08a1   listening_problem
keep folio ls health_of_child has_asthma other_chronic_illness obese listening_problem
sum health_of_child has_asthma other_chronic_illness obese listening_problem
tempfile temp_child_health
save `temp_child_health'


*******************************************************************
* Household Transfers & Expenses (bad proxy, currently not using)
*******************************************************************
* Parent transfers
use "$path/hh09dta_bx/p_tp.dta", clear
rename tp24a_2  p_g_expense_health
rename tp24b_2  p_g_other_money
rename tp24c_2  p_g_food_clothes_products
rename tp26a_2  p_r_expense_health
rename tp26b_2  p_r_other_money
rename tp26c_2  p_r_food_clothes_products
sum folio ls p_g_expense_health p_g_other_money p_g_food_clothes_products p_r_expense_health p_r_other_money p_r_food_clothes_products
* Child transfers
use "$path/hh09dta_bx/p_thi2.dta", clear
rename thi24b1_2 g_expense_health
rename thi24b2_2 g_tuition
rename thi24b4_2 g_food_clothes_products
rename thi24d1_2 r_expense_health
rename thi24d2_2 r_tuition
rename thi24d4_2 r_food_clothes_products
* Sum transfers across all reported children for each person (folio ls)
collapse (sum) g_expense_health g_tuition g_food_clothes_products r_expense_health r_tuition r_food_clothes_products, by(folio ls)
sum folio ls g_expense_health g_tuition g_food_clothes_products r_expense_health r_tuition r_food_clothes_products

*******************************************************************

********** 8888888 IMPORTANT I AM STILL FIXING THIS PROXY AS IT FALSELY SETS SOME ZEROS BUT IT'S MOSTLY CORRECT, not main variable so it's fine for now. Note: The number of observations relating to this too small. 888888 **************************


* International US Specific Remittance Proxy that includes non-cash
* Combines (a) "family member in the US" with (b) a transfer actually received from that relative 
* Only mother and father allow a clean country-of-residence + transfer match. NOTE that this is the household receives help from a member's parent living in the US, not the child's parent sends remittances. AKA could be grandparent of the child 
*******************************************************************

* Relatives in USA - yes/no (module iiib_re)
use "$path/hh09dta_b3b/iiib_re.dta", clear
keep folio ls re01
tempfile temp_re
save `temp_re'

* Relationship of relative(s) in USA (module iiib_re1)
* Restrict to spouse/father/mother/child (re03_1 1,2,3,5) to match A&C's TRN definition
* Collapsed by folio ls because a respondent can report more than one relative in the US
use "$path/hh09dta_b3b/iiib_re1.dta", clear
keep folio ls re03_1
gen close_relative_us = inlist(re03_1, 1, 2, 3, 5)
collapse (max) close_relative_us, by(folio ls)
tempfile temp_re1
save `temp_re1'

* Parent transfers + parent's country of residence (module iiib_tp)
* NOTE: mirrors the confirmed tp18pp_1/tp18pp_2 (father) coding for mother (tp18mp_1/tp18mp_2);
use "$path/hh09dta_b3b/iiib_tp.dta", clear
* We specifically look at remmittances in the united states because it has a overwhelming frequncy of 301/310, the other countries are barely any
* We code the missing appropiately. 1 is US, 0 is same location aka Mexico or a different country, . missing is categorized as don't know or genuinely missing. Note: Need to add upper in order to match across survey waves also its Estados unidos, but in 2002 it is United States
gen mother_us = .
replace mother_us = 0 if tp18mp_1 == 3
replace mother_us = 0 if tp18mp_1 == 1 & upper(tp18mp_2) != "ESTADOS UNIDOS" & !missing(tp18mp_2)
replace mother_us = 1 if tp18mp_1 == 1 & upper(tp18mp_2) == "ESTADOS UNIDOS"

gen father_us = .
replace father_us = 0 if tp18pp_1 == 3
replace father_us = 0 if tp18pp_1 == 1 & upper(tp18pp_2) != "ESTADOS UNIDOS" & !missing(tp18pp_2)
replace father_us = 1 if tp18pp_1 == 1 & upper(tp18pp_2) == "ESTADOS UNIDOS"

* TP02 = 1, parent lives in the household, need to code bc of skip logic in questionaire
* TP03 = 1, parent died more than 12 months ago, need to code bc of skip logic in questionaire
replace mother_us = 0 if tp02m == 1 | tp03m == 1
replace father_us = 0 if tp02p == 1 | tp03p == 1

* tp25 accounts for the majority of responses and that's a joint (mother & father) response. tp25m is for mothers, and tp25p is for fathers, thus we fill in receiving from mother's (recv_m) with tp25m first and tp25 if joint mother and father live together. Similarly for father
gen recv_m = tp25m
replace recv_m = tp25 if missing(recv_m) & tp22 == 1
gen recv_p = tp25p
replace recv_p = tp25 if missing(recv_p) & tp22 == 1

* remittance_us_family: 1 if either parent is confirmed in the US and sent help. 0 only if both parents are resolved as not sending. missing if unresolved and neither side is a confirmed 1 - don't silently code DK as "no". So both parents have to fail the criteria
gen remittance_us_family = .
replace remittance_us_family = 1 if (mother_us==1 & recv_m==1) | (father_us==1 & recv_p==1)
replace remittance_us_family = 0 if (mother_us==0 | (mother_us==1 & recv_m==3)) & (father_us==0 | (father_us==1 & recv_p==3))
keep folio ls remittance_us_family mother_us father_us
tempfile temp_tp_re
save `temp_tp_re'


* Combine the three pieces at the person level
use `temp_re', clear
merge 1:1 folio ls using `temp_re1', keep(master match) nogen
merge 1:1 folio ls using `temp_tp_re', keep(master match) nogen
* Genuinely a 0 because the missing are known values. They just changed to missing based on our code above. We actually know all of these and if missing it means its 0
replace close_relative_us = 0 if missing(close_relative_us)
label variable close_relative_us "Has spouse/parent/child in the US (A&C-style TRN)"
label variable remittance_us_family "Received transfer confirmed from parent living in US (A&C-style REM, parent-only)"
tempfile temp_ac_remittance
save `temp_ac_remittance'

* Collapse to household level to match the rest of the migration block
collapse (max) close_relative_us remittance_us_family, by(folio)
tempfile temp_ac_remittance_hh
save `temp_ac_remittance_hh'

*******************************************************************
* Domestic and International Remittances
* Remittances from anyone in the last 12 months outside the household in cash only. This is to at least one household member aged 15 or older 
* We use Book IIIB TP26 (parents), TH20d (siblings), THI24d (children). As the questionaire directly states these as transfers to non-residents. We account for both the survey respondents interviewed and people not interviewed via the proxy as they still are a part of the household. 
* Cash only means money for health, school tuition, or any other money. Food/clothes and time are not included
* Also we only look at the receiving end of variables as we care about individuals who have received remittances, and don't need both sides
*******************************************************************

* Parents remittances are TP26: tp26a/tp26ma/tp26pa = money for health and tp26b/tp26mb/tp26pb = any other money. We also use tp26b_1 and so forth as these are binary variables that will capture remittances if the continuous variables missed them at first. This includes mom, father, and joint. NOTE: cash_tp is all 0's and 1's and I will handle impute missing later. It is this way bc a skip pattern causes individuals who's parents died or lives in the household as skip thus they truly should be zero. 
use "$path/hh09dta_b3b/iiib_tp.dta", clear
gen cash_tp = 0
foreach v in tp26a tp26b tp26ma tp26mb tp26pa tp26pb {
    replace cash_tp = 1 if `v' > 0 & !missing(`v')
}
replace cash_tp = 1 if tp26b_1 == 1 | tp26mb_1 == 1 | tp26pb_1 == 1
* NOTE come back to this. The regressions only use the 2002 remittances, so it is left as is for now
* We account for a false zero here. TP25 aka if received help from parents has to be answered if a parent lives outside the household (tp20m or tp20p = 2). If its blank that means TP section not filled out and they are unresolved aka missing
gen miss_tp = missing(tp21) | ((tp20m == 2 | tp20p == 2) & !inlist(tp25, 1, 3) & !inlist(tp25m, 1, 3) & !inlist(tp25p, 1, 3))
keep folio ls cash_tp miss_tp
tempfile temp_cash_tp
save `temp_cash_tp'

* Sibiling remittances are TH20d. d1 = health, d2 = tuition, d3 = other money. th20c = 7 or 8 is didn't want to answer. The codebook says 7 which I believe is right, but the questionare says 8 which is a harmless addition. 
use "$path/hh09dta_b3b/iiib_th2.dta", clear
gen cash_th = (th20d1 > 0 & !missing(th20d1)) | (th20d2 > 0 & !missing(th20d2)) | (th20d3 > 0 & !missing(th20d3))
* th20c is skipped for siblings younger than 7, so a blank th20c is only a genuine skip if the sibling is under 7
gen refuse_th = inlist(th20c, 7, 8) | (missing(th20c) & !(th18_2 < 7))
collapse (max) cash_th refuse_th, by(folio ls)
tempfile temp_cash_th
save `temp_cash_th'

* Children remittances are THI24d. Same coding as siblings
use "$path/hh09dta_b3b/iiib_thi2.dta", clear
gen cash_thi = (thi24d1 > 0 & !missing(thi24d1)) | (thi24d2 > 0 & !missing(thi24d2)) | (thi24d3 > 0 & !missing(thi24d3))
* Same as siblings, thi24c is skipped for children younger than 7. A blank thi24c for a child 7 or older or age unknown is unresolved
gen refuse_thi = inlist(thi24c, 7, 8) | (missing(thi24c) & !(thi22_2 < 7))
collapse (max) cash_thi refuse_thi, by(folio ls)
tempfile temp_cash_thi
save `temp_cash_thi'

* We don't use TO aka anyone else because in 2002/2005 TO03 was only asked of people who gave help while 2009 asked everyone, so it isn't measured the same way across waves and would bias the household fixed effects

* Combine at the person level for everyone with a Book IIIB. rel = 20 is "Complete and correct", this will be another way we can check for true missing variables. Book b3b and C questionaires have this exact wording.  We rename in order not to conflict with other vars
use "$path/hh09dta_b3b/iiib_portad.dta", clear
keep folio ls rel
rename rel rel_3b
merge 1:1 folio ls using `temp_cash_tp', keep(master match) nogen
merge 1:1 folio ls using `temp_cash_th', keep(master match) nogen
merge 1:1 folio ls using `temp_cash_thi', keep(master match) nogen


* Creating new "total" variables cash3b and resolved_b to capture whethera household received remittances or whether a true missing.
foreach v in cash_tp cash_th refuse_th cash_thi refuse_thi {
    replace `v' = 0 if missing(`v')
}
egen cash_3b = rowmax(cash_tp cash_th cash_thi)
gen resolved_3b = (rel_3b == 20 & refuse_th == 0 & refuse_thi == 0 & miss_tp == 0)
keep folio ls cash_3b resolved_3b
tempfile temp_cash_3b
save `temp_cash_3b'

* We must use THE PROXY BOOK in order to account for absent household members or we will have a lot of false zeros. 
* First we do parents. Additonally we account for don't know responses. Also, made sure the variables are equal to 1 as long as they are not missing or don't know. NOTE: The naming is the same as above but the p in front is for proxy
use "$path/hh09dta_bx/p_tp.dta", clear
gen cash_ptp = 0
foreach v in tp26a tp26b tp26ma tp26mb tp26pa tp26pb {
    replace cash_ptp = 1 if (`v'_2 > 0 & !missing(`v'_2)) | (!missing(`v'_1) & `v'_1 != 8 & missing(`v'_2))
}
* DK if they don't know whether help was received TP25 = 8, or they know help was received but not what type TP26 "8. DK", stored as tp26g_1/tp26mg_1/tp26pg_1. Thus we don't know if its money or not so we code as missing. 
gen dk_ptp = (tp25 == 8 | tp25m == 8 | tp25p == 8) | !missing(tp26g_1) | !missing(tp26mg_1) | !missing(tp26pg_1)
replace dk_ptp = 1 if missing(tp21) | ((tp20m == 2 | tp20p == 2) & !inlist(tp25, 1, 3) & !inlist(tp25m, 1, 3) & !inlist(tp25p, 1, 3))
keep folio ls cash_ptp dk_ptp
tempfile temp_cash_ptp
save `temp_cash_ptp'


* We do similar thing as with siblings
use "$path/hh09dta_bx/p_th2.dta", clear
gen cash_pth = 0
foreach v in th20d1 th20d2 th20d3 {
    replace cash_pth = 1 if (`v'_2 > 0 & !missing(`v'_2)) | (!missing(`v'_1) & `v'_1 != 8 & missing(`v'_2))
}
gen dk_pth = inlist(th20c, 7, 8) | (missing(th20c) & !(th18_2 < 7))
collapse (max) cash_pth dk_pth, by(folio ls)
tempfile temp_cash_pth
save `temp_cash_pth'


* We do a similar thing with children
use "$path/hh09dta_bx/p_thi2.dta", clear
gen cash_pthi = 0
foreach v in thi24d1 thi24d2 thi24d3 {
    replace cash_pthi = 1 if (`v'_2 > 0 & !missing(`v'_2)) | (!missing(`v'_1) & `v'_1 != 8 & missing(`v'_2))
}
gen dk_pthi = inlist(thi24c, 7, 8) | (missing(thi24c) & !(thi22_2 < 7))
collapse (max) cash_pthi dk_pthi, by(folio ls)
merge 1:1 folio ls using "$path/hh09dta_bx/p_thi.dta", keepusing(thi02 thi04) nogen
replace dk_pthi = 1 if thi02 == 8 | thi04 == 8
drop thi02 thi04
tempfile temp_cash_pthi
save `temp_cash_pthi'


* Merging everything in proxy together
use "$path/hh09dta_bx/p_portad.dta", clear
keep folio ls rel
rename rel rel_px
merge 1:1 folio ls using `temp_cash_ptp', keep(master match) nogen
merge 1:1 folio ls using `temp_cash_pth', keep(master match) nogen
merge 1:1 folio ls using `temp_cash_pthi', keep(master match) nogen

gen miss_px = missing(cash_ptp)
foreach v in cash_ptp dk_ptp cash_pth dk_pth cash_pthi dk_pthi {
    replace `v' = 0 if missing(`v')
}
egen cash_px = rowmax(cash_ptp cash_pth cash_pthi)
gen resolved_px = (rel_px == 20 & dk_ptp == 0 & dk_pth == 0 & dk_pthi == 0 & miss_px == 0)
keep folio ls cash_px resolved_px
tempfile temp_cash_px
save `temp_cash_px'

* We grab the roster of the household that should contain everybody. ls01a: 1 still lives here, 4 new member, 6 returned. And we merge it with how respondents responded to receiving remittances. We use this in order to determine a true 0 rather than assign someone a false zero of no remittances. 
use "$path/hh09dta_bc/c_ls.dta", clear
* Taking care of duplicates. One is a new household member and the other is a member that left the household. We drop the person that left.
bysort folio ls: drop if _N > 1 & ls01a == 3
* We keep individuals in the household with this line
keep if inlist(ls01a, 1, 4, 6) & ls02_2 >= 15 & !missing(ls02_2)
* Merge keeps everyone by default, so anyone with a Book IIIB still counts even if they dropped out of the roster filter 
merge 1:1 folio ls using `temp_cash_3b', nogen
merge 1:1 folio ls using `temp_cash_px', keep(master match) nogen
gen cash_any = (cash_3b == 1 | cash_px == 1)
gen resolved = (resolved_3b == 1 | resolved_px == 1 | cash_any == 1)
collapse (max) cash_any (min) resolved, by(folio)

* We create our remittance variable. This is for individualks who are above the age of 15
gen remittance_hh = .
replace remittance_hh = 1 if cash_any == 1
replace remittance_hh = 0 if cash_any == 0 & resolved == 1
label variable remittance_hh "Household received remittances, proxy (indicator)" 
keep folio remittance_hh
tempfile temp_remit_hh
save `temp_remit_hh'


*******************************************************************
* Migration & Individual Remittances
*******************************************************************
* Base Income/Aid file
use "$path/hh09dta_b3a/iiia_iin.dta", clear
keep folio ls iin01d_1 iin01d_2
rename iin01d_1 ind_aid_donations_gifts
rename iin01d_2 val_aid_donations_gifts

* Merge Migration
merge 1:1 folio ls using "$path/hh09dta_b3a/iiia_mg.dta", keep(master match)
drop _merge
rename mg01p_1 country_born
rename mg04p_1 lived_upto_12

* Note this code in the preserve and restore is only adding columns and not being used right now so we don't need to worry about it for now
* Merge Migration Details collapse to person level with egen max(), revising in process
preserve
    use "$path/hh09dta_b3a/iiia_mg1.dta", clear
    * Episode-level flag: was this specific episode a move to the US?
    gen us_episode = mg30 == 1
    * Person-level: did they move to the US in ANY episode?
    sort folio ls
    by folio ls: egen moved_us = max(us_episode)
    * Person-level: documentation type reported in the US-move episode.
    gen us_doc_temp = mg31_1 if us_episode == 1
    by folio ls: egen us_doc = max(us_doc_temp)
    sort folio ls
    by folio ls: keep if _n == _N
    keep folio ls moved_us us_doc
    tempfile temp_mg1
    save `temp_mg1'
restore
merge 1:1 folio ls using `temp_mg1', keep(master match) nogen

* This equivalent coding to using val_aid_donations_gits > 1, this is currently encompassing a lot more than remittances and its possible it doesn't encompass remittances at all 
gen remittance_binary = .
replace remittance_binary = 1 if ind_aid_donations_gifts == 1
replace remittance_binary = 0 if ind_aid_donations_gifts == 3
label variable remittance_binary "Household Received Aid/Donations (1=Yes, 0=No)"
drop ind_aid_donations_gifts
tempfile temp_migration

* We collapse to household level because these variables are for adults, thus to match with the kids biomarkers we need to collapse. Thus, we do so and do a many to one merge
collapse (max) remittance_binary val_aid_donations_gifts country_born lived_upto_12 moved_us us_doc, by(folio)
tempfile temp_migration_hh
save `temp_migration_hh'

*******************************************************************
* We need linking key across waves
*******************************************************************

* We specifically take the pid_link in this one because it has the most nonmissing pid_link rows, on second thought I don't think it needs to be this one, but it doesn't hurt to do so
use "$path/hh09dta_bc/c_ls.dta", clear
* Taking care of duplicates. One is a new household member and the other is a member that left the household. We drop the person that left.
bysort folio ls: drop if _N > 1 & ls01a == 3
keep folio ls pid_link

tempfile temp_pidlink
save `temp_pidlink'

*******************************************************************
* Demographics & Control Variables
*******************************************************************

* We keep only the state code (ent) from the Book IIIA cover, collapsed to the household, because the state fill-in block further down uses it. Same file, rows and collapse as before, so ent is unchanged
use "$path/hh09dta_b3a/iiia_portad.dta", clear
keep folio ent
collapse (max) ent, by(folio)
tempfile temp_controls_hh
save `temp_controls_hh'

*******************************************************************
* Child-Level Controls: Vaccination, Delivery & Breastfeeding
*******************************************************************

* Child Vaccination (Book 5, v_vac), this is a clean individual child level control
* NOTE THIS CONTROL IS RATHER CUMBERSOME
use "$path/hh09dta_b5/v_vac.dta", clear
* Does the child have a current vaccination card, a previous card, dates copied off the card, and then a recall quesiton asking the parents if their child was vaccinated
keep folio ls vac02 vac02_a vac03a_1-vac03x_1 vac04a vac04b vac04c vac04d vac04e vac04f vac04g vac04h
rename vac02 has_vaccination_card
* IMPORTANT NOTE: Skip pattern, VAC04 is asked ONLY when no card is seen. If a card is seen VAC02 = 1, or 2005-08 card VAC02_a = 1 the interviewer copies the dates into VAC03 and skips VAC04.

* This generates a categorical variable that denotes if 1 = you had a card, 2 = a previous card, or 3 = a recall.
gen byte card_src = 1 if has_vaccination_card == 1
replace card_src = 2 if vac02_a == 1 & missing(card_src)
replace card_src = 3 if inlist(has_vaccination_card, 2, 3)

* Writing a loop to check if there is a date corresponding with if the child is vaccinated. So r_ = 1 if the card shows the dose was given a real date, or "44" = given but date not recorded, and 0 if the box is blank or "45". 45 means not vaccinateed because too young
foreach v of varlist vac03a_1-vac03x_1 {
    gen byte r_`v' = (trim(`v') != "" & substr(trim(`v'), 1, 2) != "45")
}

* Combining doses into vaccines, so if a child got any of these doses and their 2009 card was seen (card_src = 1) then it is 1. Pentavalent acellular = DPaT + IPV (polio) + Hib
gen byte v_penta = max(r_vac03e_1, r_vac03f_1, r_vac03g_1, r_vac03h_1) if card_src == 1
gen byte v_bcg = r_vac03a_1 if card_src == 1
gen byte v_polio = max(v_penta, r_vac03t_1) if card_src == 1
gen byte v_dpt = max(v_penta, r_vac03i_1) if card_src == 1
gen byte v_mmr = max(r_vac03r_1, r_vac03s_1) if card_src == 1
gen byte v_measles = max(v_mmr, r_vac03v_1) if card_src == 1
gen byte v_hepb = max(r_vac03b_1, r_vac03c_1, r_vac03d_1, r_vac03w_1, r_vac03x_1) if card_src == 1


* Same idea but now we replace them if their previous card shows a dosage.
replace v_penta = max(r_vac03i_2, r_vac03j_2, r_vac03k_2) if card_src == 2
replace v_bcg = r_vac03a_2 if card_src == 2
replace v_polio = max(r_vac03b_2, r_vac03c_2, r_vac03d_2, r_vac03e_2, r_vac03f_2, r_vac03g_2, r_vac03h_2) if card_src == 2
replace v_dpt = max(v_penta, r_vac03l_2, r_vac03m_2, r_vac03n_2, r_vac03o_2, r_vac03p_2) if card_src == 2
replace v_mmr = max(r_vac03r_2, r_vac03s_2) if card_src == 2
replace v_measles = max(v_mmr, r_vac03q_2) if card_src == 2
replace v_hepb = max(v_penta, r_vac03t_2) if card_src == 2


* For recall (card_src = 3), we use a nesting condition. Basically if parent answered yes then 1, no(3) then 0, and don't know or missing is set to missing. Same logic as 2005 cleaning script
replace v_bcg = cond(vac04a == 1, 1, cond(vac04a == 3, 0, .)) if card_src == 3
replace v_polio = cond(vac04b == 1, 1, cond(vac04b == 3, 0, .)) if card_src == 3
replace v_dpt = cond(vac04c == 1, 1, cond(vac04c == 3, 0, .)) if card_src == 3
replace v_penta = cond(vac04d == 1, 1, cond(vac04d == 3, 0, .)) if card_src == 3
replace v_mmr = cond(vac04e == 1, 1, cond(vac04e == 3, 0, .)) if card_src == 3
replace v_measles = cond(vac04f == 1 | vac04e == 1, 1, cond(vac04f == 3 & vac04e == 3, 0, .)) if card_src == 3
replace v_hepb = cond(vac04g == 1, 1, cond(vac04g == 3, 0, .)) if card_src == 3


* We now use egen to calulate the number of vaccines taken by a child. Fully vaccinated is 0 if any of the 7 is known NOT received. 1 if all 7 received. Missing otherwise
egen byte vac_min7 = rowmin(v_bcg v_polio v_dpt v_penta v_mmr v_measles v_hepb)
egen byte vac_nmiss7 = rowmiss(v_bcg v_polio v_dpt v_penta v_mmr v_measles v_hepb)
gen fully_vaccinated = 1 if vac_min7 == 1 & vac_nmiss7 == 0
replace fully_vaccinated = 0 if vac_min7 == 0
label variable fully_vaccinated "Fully vaccinated (indicator)" 
drop r_vac03* vac03* vac_nmiss7 vac_min7 v_* card_src
tempfile temp_vaccination
save `temp_vaccination'



* The 2009 roster has 8 line numbers shared by a person who left and a current person, so this drops the person who left, the same fix the script uses at the top 
use "$path/hh09dta_bc/c_ls.dta", clear
bysort folio ls: drop if _N > 1 & ls01a == 3
tempfile temp_roster09
save `temp_roster09'


* NOTE: THIS IS A HUGE DIFFERENCE TO 2002 bc no child line number recorded thus we have to do it ourselves. 
* NOTE: Unlike 2002, the 2005 pregnancy records don't say which child each pregnancy belongs to. So we link them using the mother's line number plus the baby's sex and birth month and year (or birth order when no dates exist), and drop any link that isn't unique. Children also in 2002 then get their exact 2002 values, so the matching mainly fills in children born after 2002.

* This codechunk is grabbing the child's characteristics and only keeping mothers that are in the household and not dead. Also we make both child_ls and mother_ls. Note child_ls is not just for children it's for all individuals but we name it child_ls because we are going to use it for matching
use `temp_roster09', clear
keep folio ls ls07 ls03_21 ls03_22 ls04
gen mother_ls = string(ls07, "%02.0f") if ls07 < 51
drop if mother_ls == ""
rename ls child_ls
tempfile temp_kids_mother
save `temp_kids_mother'
* Here we are calculating the full histories. The first codechunk is calculating household, mother's line, and number of pregnancies belonging to a mother
use "$path/hh09dta_b4/iv_he1.dta", clear
bysort folio ls: egen n_preg = max(secuencia)
preserve
    keep folio ls n_preg
    duplicates drop
    rename ls mother_ls
    tempfile temp_he1_women
    save `temp_he1_women'
restore
keep if he10b == 2
sort folio ls secuencia
by folio ls: gen byte birth_order = _n
gen byte oldest_child = (birth_order == 1)
gen byte delivered_home = inlist(he15_1, 10, 11, 12) if !missing(he15_1)
keep folio ls secuencia birth_order oldest_child delivered_home he11 he12_22 he12_23
rename ls mother_ls
* Link to the roster child with the same mother, birth year, birth month and sex. A pregnancy matching more than one child, or a child matching more than one pregnancy, is dropped
* NOTE: THIS FIRST MERGE GIVES US birth order, oldest child and delivered at home, the variables that live in the full history.
joinby folio mother_ls using `temp_kids_mother'
keep if he12_23 == ls03_22 & he12_22 == ls03_21 & he11 == ls04 & !missing(he12_23, he12_22, he11)
bysort folio mother_ls secuencia: keep if _N == 1
bysort folio child_ls: keep if _N == 1
keep folio mother_ls secuencia child_ls birth_order oldest_child delivered_home
tempfile temp_he1_link
save `temp_he1_link'

* Baby birthweight
use "$path/hh09dta_b4/iv_he2.dta", clear
* Only keeping babies that were born alive 
keep if he20b == 2
* 2005 difference: birth weight he46_2 is already in kg, while 2002 was in grams. Codebook is wrong for 2005 it makes it seem like 2005 is in grams but the data shows that it is not
gen birth_weight_kg = he46_2 if he46_1 == 1
* Implausible weights are set to missing, we air on the side of caution and have a large upper bound of 8 kg or more
replace birth_weight_kg = . if birth_weight_kg >= 8

* Weeks for predominantly breastfed. Meaning milk, tea, or water. Typically this could be shorter than just total breastfed which no variable exists here to capture the duration
gen predominant_breastfed = he48_2 / 7 if he48_1 == 1
replace predominant_breastfed = he48_3 if he48_1 == 2
replace predominant_breastfed = he48_4 * 30.4375 / 7 if he48_1 == 3
replace predominant_breastfed = 0 if he47 == 3
gen byte still_breastfed = (he48_1 == 4) if !missing(he47)

* Iron during pregnancy
gen iron_pregnancy = .
replace iron_pregnancy = 1 if he26a == 1
replace iron_pregnancy = 0 if he26a == 3

* Premature birth
gen premature = .
replace premature = 1 if he30 == 1
replace premature = 0 if he30 == 3

* We build a second delivered at home for women with no full history because they don't have a he15_1
gen byte delivered_home_he2 = inlist(he31_1, 10, 11, 12) if inrange(he31_1, 1, 13) | he31_1 == 15
rename ls mother_ls

* Keeping the relevant variables and then merging in order to label if the mother is matched or not. Which 1 is matched and 0 is not
keep folio mother_ls secuencia birth_weight_kg predominant_breastfed still_breastfed iron_pregnancy premature delivered_home_he2
merge m:1 folio mother_ls using `temp_he1_women', keep(master match)
gen byte has_he1 = (_merge == 3)
drop _merge
tempfile temp_he2_vars
save `temp_he2_vars'



*** REVISING THIS CODE BELOW UNTIL THE END OF THE SECTION, BUT NOT USING IN DID RN


* NOTE: THIS SECOND MERGE GIVES US birth weight, breastfeeding, iron and premature, aren't in the full history. They're in a different file, the recent-pregnancy file
* Step A
* We keep only women who have a full history otherwise they go to Step B
keep if has_he1 == 1
* Converts each recent record's number into the full-history number.
gen secuencia_he1 = n_preg - secuencia + 1
drop secuencia
rename secuencia_he1 secuencia
* Swaps in the converted number under the same name as in temp_he1_link, so the two files can be merged on it.
merge 1:1 folio mother_ls secuencia using `temp_he1_link', keepusing(child_ls) keep(match) nogen
tempfile temp_he2_a
save `temp_he2_a'

* Step B
* Handle mothers who skipped the full history
use `temp_kids_mother', clear
* We only keep babies after 2005
keep if ls03_22 >= 2005 & !missing(ls03_22, ls03_21)
* gsort with minus signs: sorts each mother's children by birth year, then month, newest first, the minus means descending. Then we count the kids rank aka k and the number of kids
gsort folio mother_ls -ls03_22 -ls03_21
by folio mother_ls: gen k = _n
by folio mother_ls: gen n_kids = _N
keep folio mother_ls k n_kids child_ls
tempfile temp_kids_order
save `temp_kids_order'
use `temp_he2_vars', clear
keep if has_he1 == 0
sort folio mother_ls secuencia
* Counts how many live births since 2005 she reported 
by folio mother_ls: gen k = _n
by folio mother_ls: gen n_births = _N
* Pairs birth k with child k for the same mother: newest birth to youngest child, second-newest to second-youngest.
merge 1:1 folio mother_ls k using `temp_kids_order', keep(match) nogen
* Only keep a mother's links if she reported exactly as many births as she has children born since 2005 on the roster 
keep if n_births == n_kids
append using `temp_he2_a'
* A child linked twice is ambiguous, so both links are dropped
bysort folio child_ls: drop if _N > 1
keep folio child_ls birth_weight_kg predominant_breastfed still_breastfed iron_pregnancy premature delivered_home_he2
tempfile temp_he2_link
save `temp_he2_link'

* Combines the two link tables into one file with one row per child and all the Book IV variables, ready to merge onto the 2005 children
use `temp_he1_link', clear
keep folio child_ls birth_order oldest_child delivered_home
merge 1:1 folio child_ls using `temp_he2_link', nogen
replace delivered_home = delivered_home_he2 if missing(delivered_home)
drop delivered_home_he2
rename child_ls ls
label variable birth_order "Birth order" 
label variable oldest_child "Oldest child (indicator)" 
label variable birth_weight_kg "Birth weight (kilograms)"
label variable delivered_home "Delivered at home (indicator)" 
label variable predominant_breastfed "Weeks predominantly breastfed"
label variable still_breastfed "Still being breastfed at interview"
label variable iron_pregnancy "Iron during pregnancy (indicator)"
label variable premature "Premature, mother reported (indicator)"
tempfile temp_bookiv_child
save `temp_bookiv_child'


*******************************************************************
* Parent Controls: at the child level 
*******************************************************************

* Mother's and father's height. We use measured height, not self-reported height to match the child's height. Overarching idea is that book s has parents heights but it doesn't tell us who the childs parent is, thus we merge it to book c which does
* IMPORTANT NOTE: Father's height is a bit sketchy bc a lot of them are missing. And this could be due to migration. Thus, it might be better to not use it in the regression. Could omit both and see the results as well.
* The difference between wave 2002 is that it's sa07_21 not sa07_2 and it's ls not ls00
use "$path/hh09dta_bs/s_sa.dta", clear
keep folio ls sa03 sa07_21
* The 2 duplicated Book S people are split records of the same person (see the top of the script), so we keep the row that has the height
bysort folio ls (sa07_21): keep if _n == 1
* We must make sure that the book S person is the same person by checking if their age is within 5 years, we also keep gender for later
merge 1:1 folio ls using `temp_roster09', keepusing(ls02_2 ls04) keep(master match) nogen
replace sa07_21 = . if abs(sa03 - ls02_2) > 5 | missing(sa03, ls02_2)
* We set heights that are not possible to missing, only a handful
replace sa07_21 = . if !inrange(sa07_21, 120, 210)
keep folio ls sa07_21 ls04
tempfile temp_heights
save `temp_heights'
* We replace observations who correspond to the wrong sex
replace sa07_21 = . if ls04 != 3
drop ls04
rename (ls sa07_21) (mother_ls mother_height_cm)
tempfile temp_mother_height
save `temp_mother_height'
use `temp_heights', clear
replace sa07_21 = . if ls04 != 1
drop ls04
rename (ls sa07_21) (father_ls father_height_cm)
tempfile temp_father_height
save `temp_father_height'
* We merge in the parents heights to their children. Also not that we use less than 51 because that means the parent is dead or not living in the household. This is put here just bc this is missing data. We understand that it doesn't matter if a parent is dead or alive to contribute to child height. Note ls is a string for 2005 wave but is a number for 2002
use `temp_roster09', clear
keep folio ls ls06 ls07
gen mother_ls = string(ls07, "%02.0f") if ls07 < 51
gen father_ls = string(ls06, "%02.0f") if ls06 < 51
merge m:1 folio mother_ls using `temp_mother_height', keep(master match) nogen
merge m:1 folio father_ls using `temp_father_height', keep(master match) nogen
keep folio ls mother_height_cm father_height_cm
label variable mother_height_cm "Mother's height (centimeters)" 
label variable father_height_cm "Father's height (centimeters)" 
tempfile temp_parent_height
save `temp_parent_height'

****************************************************************************
* TO DO NOTE: Perhaps add mother's and father's WEIGHT, currently not done
****************************************************************************


* Mother's and father's age at the child's birth, in years. Computed as the parent's roster age minus the child's roster age, since MxFLS 2005 has no direct question for this. Roster ages are in whole years, so this can be off by up to 1 year
use `temp_roster09', clear
keep folio ls ls02_1 ls02_2 ls04
* 2009 difference: ls02_1 is only 1 = age given in years or 8 = DK, there are no months or days units like 2002.
gen age_years = ls02_2 if ls02_1 == 1
keep folio ls age_years ls04
tempfile temp_ages
save `temp_ages'
* The gender must match across files. Similar to the code above for height
replace age_years = . if ls04 != 3
drop ls04
rename (ls age_years) (mother_ls mother_age)
tempfile temp_mother_age
save `temp_mother_age'
use `temp_ages', clear
replace age_years = . if ls04 != 1
drop ls04
rename (ls age_years) (father_ls father_age)
tempfile temp_father_age
save `temp_father_age'
* Link each child to their mother and father. 51 = parent dead or not living in the household, so their age is not known and it stays missing, not 0. Unlike 2002 wave, here the ls is a character 
use `temp_roster09', clear
keep folio ls ls06 ls07
merge 1:1 folio ls using `temp_ages', keepusing(age_years) nogen
gen mother_ls = string(ls07, "%02.0f") if ls07 < 51
gen father_ls = string(ls06, "%02.0f") if ls06 < 51
merge m:1 folio mother_ls using `temp_mother_age', keep(master match) nogen
merge m:1 folio father_ls using `temp_father_age', keep(master match) nogen
gen mother_age_at_birth = mother_age - age_years
gen father_age_at_birth = father_age - age_years
* Implausible ages at birth are set to missing, same bounds as 2002
replace mother_age_at_birth = . if !inrange(mother_age_at_birth, 10, 60)
replace father_age_at_birth = . if !inrange(father_age_at_birth, 10, 80)
keep folio ls mother_age_at_birth father_age_at_birth
label variable mother_age_at_birth "Mother's age at birth (years)" 
label variable father_age_at_birth "Father's age at birth (years)"
tempfile temp_parent_age
save `temp_parent_age'


* Mother's and father's employment last week
use "$path/hh09dta_b3a/iiia_tb.dta", clear
keep folio ls tb02_1 tb03 tb04 tb05
* We use proxy to account for individuals that were not present at interview. This saves over 1,200 fathers. We only keep one row per individual
append using "$path/hh09dta_bx/p_tb.dta", keep(folio ls tb02_1 tb03 tb04 tb05) gen(proxy)
gen employed = 0 if tb05 == 3
replace employed = 1 if tb02_1 == 1 | tb03 == 1 | tb04 == 1 | tb05 == 1
* We only keep one row per individual, we sort so we have priority of Book IIIA with answer, Proxy with answer, Book IIIA with no answer, and then proxy with no answer
gen byte noans = missing(employed)
sort folio ls noans proxy
by folio ls: keep if _n == 1
drop noans
* We merge to bring in the sex checks
merge 1:1 folio ls using `temp_roster09', keepusing(ls04) keep(master match) nogen
keep folio ls employed ls04
tempfile temp_employed
save `temp_employed'
replace employed = . if ls04 != 3
drop ls04
rename (ls employed) (mother_ls mother_employed)
tempfile temp_mother_employed
save `temp_mother_employed'
use `temp_employed', clear
replace employed = . if ls04 != 1
drop ls04
rename (ls employed) (father_ls father_employed)
tempfile temp_father_employed
save `temp_father_employed'
* Link each child to their mother and father. 51 = parent dead or not living in the household, so employment is not known and stays missing, not 0. A parent under 15 or with no Book IIIA or Proxy answer also stays missing. Unlike 2002 wave, here the ls is a character 
use `temp_roster09', clear
keep folio ls ls06 ls07
gen mother_ls = string(ls07, "%02.0f") if ls07 < 51
gen father_ls = string(ls06, "%02.0f") if ls06 < 51
merge m:1 folio mother_ls using `temp_mother_employed', keep(master match) nogen
merge m:1 folio father_ls using `temp_father_employed', keep(master match) nogen
keep folio ls mother_employed father_employed
label variable mother_employed "Mother employed (indicator)" 
label variable father_employed "Father employed (indicator)"
tempfile temp_parent_employed
save `temp_parent_employed'


* Mother's and father's years of education completed in years
* We add the proxy book for household members not present. Also we rename these variables bc they are the exact same in the b3 not proxy.
use "$path/hh09dta_bx/p_ed.dta", clear
keep folio ls ed04 ed05 ed06_1 ed07
rename (ed04 ed05 ed06_1 ed07) (ed05 ed06 ed07_1 ed08)
tempfile temp_proxy_ed
save `temp_proxy_ed'
use "$path/hh09dta_b3a/iiia_ed.dta", clear
keep folio ls ed05 ed06 ed07_1 ed08
* Same as 2002 and 2005, we drop duplicates only after edu_years is built, keeping the row that has an answer. 553 people are in both books in 2009 and for 438 of them the Book IIIA section is blank
append using `temp_proxy_ed', gen(proxy)
* We code ED05 = 3 (never attended school) and ED06 = 1 or 2 (no instruction, preschool) as 0 years
gen edu_years = 0 if ed05 == 3 | inlist(ed06, 1, 2)
* First line is for elementary school and is just the value associate ed with ed07_1
* Second line is elementary school which is 6 years in mexico plus the additional years of middle school that goes up until 3. If it is above 3 then it is set to missing. 
* Third line works the same but for highschool. 
replace edu_years = ed07_1 if ed06 == 3 & inrange(ed07_1, 0, 6)
replace edu_years = 6 + ed07_1 if ed06 == 4 & inrange(ed07_1, 0, 3)
replace edu_years = 9 + ed07_1 if ed06 == 6 & inrange(ed07_1, 0, 3)
* Here the open highschool, open secondary, normal basic (which means high school in the questionaire), college, and graduate are not asked how many years completed rather we are just given a metric for complted or not aka we skip ED07 and only ask ED08, whether they got the degree.
* Also note we take the lower bound, if not completed we basically say they just finished the previous degree and started the new degree (e.g. highschool not completed means they get coded a completed middle school level education). Similar strategy with college.
replace edu_years = 9 if ed06 == 5 & ed08 == 1
replace edu_years = 6 if ed06 == 5 & inlist(ed08, 2, 3)
replace edu_years = 12 if inlist(ed06, 7, 8) & ed08 == 1
replace edu_years = 9 if inlist(ed06, 7, 8) & inlist(ed08, 2, 3)
replace edu_years = 16 if ed06 == 9 & ed08 == 1
replace edu_years = 12 if ed06 == 9 & inlist(ed08, 2, 3)
* We take 18 years if they had graduate school as this is the most common graduate program. It's possible to be 1 year, but I think 2 years is a safe choice still. Lower bound for PhD
replace edu_years = 18 if ed06 == 10 & ed08 == 1
replace edu_years = 16 if ed06 == 10 & inlist(ed08, 2, 3)
* Drop duplicates only now that edu_years is built, keeping the row that has an answer
gen byte noans = missing(edu_years)
sort folio ls noans proxy
by folio ls: keep if _n == 1
drop noans
merge 1:1 folio ls using `temp_roster09', keepusing(ls04) keep(master match) nogen
keep folio ls edu_years ls04
tempfile temp_edu_years
save `temp_edu_years'
* Making sure gender match for mother and father. Same as the code above. 
replace edu_years = . if ls04 != 3
drop ls04
rename (ls edu_years) (mother_ls mother_edu_years)
tempfile temp_mother_edu
save `temp_mother_edu'
use `temp_edu_years', clear
replace edu_years = . if ls04 != 1
drop ls04
rename (ls edu_years) (father_ls father_edu_years)
tempfile temp_father_edu
save `temp_father_edu'
* Link each child to their mother and father. 51 = parent dead or not living in the household, so education is not known and stays missing, not 0. Unlike 2002 wave, here the ls is a character 
use `temp_roster09', clear
keep folio ls ls06 ls07
gen mother_ls = string(ls07, "%02.0f") if ls07 < 51
gen father_ls = string(ls06, "%02.0f") if ls06 < 51
merge m:1 folio mother_ls using `temp_mother_edu', keep(master match) nogen
merge m:1 folio father_ls using `temp_father_edu', keep(master match) nogen
keep folio ls mother_edu_years father_edu_years
label variable mother_edu_years "Mother's years of schooling" 
label variable father_edu_years "Father's years of schooling" 
tempfile temp_parent_edu
save `temp_parent_edu'



*******************************************************************
* Household-Level Controls: Household Head & Dwelling
*******************************************************************


* Dwelling characteristics: Electricity and household floor are proxies of wealth
use "$path/hh09dta_bc/c_cvo.dta", clear
keep folio cvo04 cvo05_1 cvo07_1 cvo03a cvo03b cvo03c cvo03d cvo03e cvo03f

* Electricity (binary)
gen dwell_electricity = (cvo04 == 1) if !missing(cvo04)
label variable dwell_electricity "Dwelling has electricity (indicator)" 

* Good floor quality
gen dwell_floor_good = inlist(cvo05_1, 1, 2) if inrange(cvo05_1, 1, 4)
label variable dwell_floor_good "Good floor quality (indicator)"

* Good roof quality = 1 if beam/polyurethane, concrete/brick, slate/tile, metal/fiberglass/plastic sheets, or palm/shingle/wood (codes 1, 2, 3, 6, 7). 0 if asbestos, reed-grass/bamboo, cardboard, waste material or other (codes 4, 5, 8, 9, 10) 
gen dwell_roof_good = inlist(cvo07_1, 1, 2, 3, 6, 7) if inrange(cvo07_1, 1, 10)
label variable dwell_roof_good "Good roof quality (indicator)" 

* Good ventilation. We grab all the HH dwelling 6 options so if its present we can code it as 0 and if all missing we can code it as missing. ==4 is the one about ventilation. And the 3rd line of code corrects for a contradiction btw none of the above and dwelling has enough ventilation
egen cvo03_n = rownonmiss(cvo03a cvo03b cvo03c cvo03d cvo03e cvo03f)
gen dwell_ventilation_good = (cvo03d == 4) if cvo03_n > 0
replace dwell_ventilation_good = . if cvo03d == 4 & cvo03f == 6
label variable dwell_ventilation_good "Enough ventilation (indicator)" 

*** Okay we have about 8% or 850 missing floor values, doesn't seem to be a skip pattern or anything, I was going to match it with previous survey waves but decided not to as this might assume/understate if they have upgraded, thus if its balanced across the remittances it should not bias WE CAN CHECK THIS BALANCE I GUESS. 
keep folio dwell_*
tempfile temp_dwell_cvo
save `temp_dwell_cvo'

* Rooms per person is rooms used for sleeping / household size. 
use `temp_roster09', clear

* Roster members who currently live in the household. We filter to only ppl living in the household, new and old household members (ls01a = 1 still lives here, 4 new member, 6 returned). ls09 is not in this wave unlike 2002. Thus, we count rows with a 1 instead of counting ls. Also slightly different from wave 2005 as it has 6 for returned
keep if inlist(ls01a, 1, 4, 6)
gen byte one = 1
collapse (sum) hh_size = one, by(folio)
* CV07 is asked of every household. 0 rooms is a real answer as 1 household, they sleep in the kitchen. 14 households report exactly 10 rooms although they have only 1 to 6 members and only 1 household reports 9, so 10 looks like a code or typing error and we set it to missing
* We cap households with rooms for sleeping under 10, this is reasonable gets rid of these huge outliers
merge 1:1 folio using "$path/hh09dta_bc/c_cv.dta", keepusing(cv07) keep(master match) nogen
replace cv07 = . if cv07 >= 10 & !missing(cv07)
gen rooms_per_person = cv07 / hh_size
label variable rooms_per_person "Sleeping rooms per household member" 
keep folio rooms_per_person
tempfile temp_rooms
save `temp_rooms'

* Monthly per capita household expenditure. Total monthly household expenditure / household size, in pesos per person per month
* Note we don't include gifts given away or school spending on non fimily members
use "$path/hh09dta_b1/i_cs.dta", clear
merge 1:1 folio using "$path/hh09dta_b1/i_cs1.dta", nogen
foreach f of varlist cs02?_?1 cs04?_*1 cs16?_1 cs18_1 cs22?_1 cs24?_1 cs27?_1 cs29_1 cs34a_11 cs34a_21 cs35a_11 cs35a_21 cs36a_11 cs36a_21 {
    local v = substr("`f'", 1, strlen("`f'") - 1) + "2"
    replace `v' = 0 if inlist(`v', 99, 999, 9999, 99999)
    gen byte dk_tmp = (`f' == 8 | inlist(`v', 98, 998, 9998, 99998))
    replace `v' = . if dk_tmp == 1
    egen med_tmp = median(`v')
    replace `v' = med_tmp if dk_tmp == 1
    drop med_tmp dk_tmp
}
* We set cs15 to missing if the household said no they didn't buy that good aka cs08 isn't 1
foreach L in a b c d e f g h i j {
    replace cs15`L' = . if cs08`L' != 1
    replace cs15`L' = . if inlist(cs15`L', 98, 99, 998, 999, 9998, 9999, 99998, 99999)
    egen med_tmp = median(cs15`L')
    replace cs15`L' = med_tmp if cs08`L' == 1 & missing(cs15`L')
    drop med_tmp
}

* We add up each row, then convert every unit to months
egen food_week = rowtotal(cs02?_?2 cs15? cs04?_*2)
egen month_items = rowtotal(cs16?_2 cs18_2)
egen quarter_items = rowtotal(cs22?_2 cs24?_2)
egen year_items = rowtotal(cs27?_2 cs29_2)
egen school_year = rowtotal(cs34a_12 cs34a_22 cs35a_12 cs35a_22)
egen school_transport = rowtotal(cs36a_12 cs36a_22)
* School transport is converted to years and then months, same with food week. Also the school_year is just divided by 12 for months not by 10.5 as we spread out the cost (also we already double checked why we are using this 12 instead of 10.5)
gen hh_expend_month = food_week * 52 / 12 + month_items + quarter_items / 3 + year_items / 12 + school_year / 12 + school_transport * 10.5 / 12
* Rowtotal counts blanks as 0, which is right for items that were skipped, but 14 households did not finish the section - a filter question is blank, so their total stays missing not 0
replace hh_expend_month = . if missing(cs03, cs17, cs23, cs28, cs32, cs08a, cs08b, cs08c, cs08d, cs08e, cs08f, cs08g, cs08h, cs08i, cs08j)
keep folio hh_expend_month
tempfile temp_expend
save `temp_expend'
* Household size to create the monthly expenditure for a household. This uses ls01a instead of ls09 like wave 2002. Also adds a 6 returned to household unlike wave 2005
use `temp_roster09', clear
keep if inlist(ls01a, 1, 4, 6)
gen byte one = 1
collapse (sum) hh_size = one, by(folio)
merge 1:1 folio using `temp_expend', keep(match) nogen
gen pc_expend_month = hh_expend_month / hh_size
label variable pc_expend_month "Monthly per capita household expenditure in pesos"
* Log of per capita expenditure for the regressions since the level is very skewed
gen ln_pc_expend_month = ln(pc_expend_month)
label variable ln_pc_expend_month "Ln (monthly expenditure per capita, pesos), proxy" 
keep folio pc_expend_month ln_pc_expend_month
tempfile temp_pc_expend
save `temp_pc_expend'

***************
** TESTING ****
***************
* Maternal/Household Child Mortality (Book 4, Reproductive Health)
* NOTE: iv_res is a mother-level LIFETIME cumulative measure (every child she's ever borne, as of the 2009 interview) -- not a birth-year-specific outcome like. We collapse it to household level
*use "$path/hh09dta_b4/iv_res.dta", clear
*keep folio ls res03 res04 res07 res08 res10 res11 res13
*rename res10 sons_died
*rename res11 daughters_died
*rename res13 stillbirths
* children_ever_born_alive = living with her (res03/04) + alive elsewhere (res07/08)
*                            + died after being born alive (sons_died/daughters_died)
*gen children_ever_born_alive = res03 + res04 + res07 + res08 + sons_died + daughters_died
*gen children_died = sons_died + daughters_died
*tempfile temp_mortality
*save `temp_mortality'
* Sum numerator and denominator ACROSS all mothers in the household before
* dividing -- summing first avoids a "ratio of ratios" mistake if a household
* has more than one woman with a birth history.
*collapse (sum) children_ever_born_alive children_died stillbirths, by(folio)
*gen child_mortality_rate = 100 * children_died / children_ever_born_alive if children_ever_born_alive > 0
*label variable child_mortality_rate "HH lifetime child deaths per 100 live births (mothers' birth histories, cumulative to 2009)"
*tempfile temp_mortality_hh
*save `temp_mortality_hh'
***************
** TESTING ****
***************

*******************************************************************
* Merging all datasets
*******************************************************************
use `temp_biomarkers', clear
merge 1:1 folio ls using `temp_vaccination', keep(master match) gen(_merge_vaccination)
merge m:1 folio using `temp_controls_hh', keep(master match) gen(_merge_controls)
merge m:1 folio using `temp_migration_hh', keep(master match) gen(_merge_migration)
merge m:1 folio using `temp_ac_remittance_hh', keep(master match) gen(_merge_ac_remittance)
merge m:1 folio using `temp_remit_hh', keep(master match) gen(_merge_remit_hh)
merge 1:1 folio ls using `temp_pidlink', keep(master match) gen(_merge_pidlink)
merge m:1 folio using `temp_dwell_cvo', keep(master match) gen(_merge_dwell)
merge 1:1 folio ls using `temp_bookiv_child', keep(master match) gen(_merge_bookiv)
* Same as 2002: a child still being breastfed at the interview has been predominantly breastfed for at least their age, so we use their age in weeks
replace predominant_breastfed = age_days / 7 if still_breastfed == 1 & age_days >= 0 & !missing(age_days)
* Predominant breastfeeding over 5 years (260.9 weeks) is not plausible. This is a conservative cutoff, only about 15-20 children per wave are above it
replace predominant_breastfed = . if predominant_breastfed > 5 * 365.25 / 7 & !missing(predominant_breastfed)
merge 1:1 folio ls using `temp_parent_height', keep(master match) gen(_merge_parent_height)
merge 1:1 folio ls using `temp_parent_age', keep(master match) gen(_merge_parent_age)
merge 1:1 folio ls using `temp_parent_employed', keep(master match) gen(_merge_parent_employed)
merge 1:1 folio ls using `temp_parent_edu', keep(master match) gen(_merge_parent_edu)
merge m:1 folio using `temp_rooms', keep(master match) gen(_merge_rooms)
merge m:1 folio using `temp_pc_expend', keep(master match) gen(_merge_pc_expend)

* Ent only came from the Book IIIA cover, iiia_portad, in the controls collapse, so a household where no adult has a Book IIIA got a missing state and was dropped from the rainfall merge. We fill this in with this book below. This matches the other waves
preserve
    use "$path/hh09dta_bc/c_portad.dta", clear
    keep folio ent
    rename ent ent_c
    tempfile temp_ent_c
    save `temp_ent_c'
restore
merge m:1 folio using `temp_ent_c', keep(master match) nogen
replace ent = ent_c if missing(ent)
drop ent_c

*******************************************************************
* Municipality
*******************************************************************

* We keep the Book C state aka edo together with mpio so each municipality sits inside its own state, and build the 5 digit INEGI key cvegeo = 2 digit state + 3 digit municipality.
* NOTE: 21 Book C cover rows have no state or municipality, all with book interview result rel = 21. None of these households has a measured child, so we leave them missing 
preserve
    use "$path/hh09dta_bc/c_portad.dta", clear
    keep folio ent mpio loc
	* It was checked for the majority that 35 households coded 15-011 in 2009 come from households that were in 15-110 in 2002 and 2005. Contrary to the mg08a survey that says these households didn't move. That accounts for most of them but a couple are assumed, but this is a legit assumption. So 011 is a typo for 110, digits reversed.
    replace mpio = 110 if ent == 15 & mpio == 11
    * We checked this. This is likely coded incorrectly because same locality, same January 2010 field trip, and the MG08a answers. Matches the other survey waves
    replace mpio = 102 if ent == 16 & mpio == 10 & loc == 15
    gen str5 cvegeo = string(ent, "%02.0f") + string(mpio, "%03.0f") if !missing(ent, mpio)
    label variable mpio "Municipality of residence, INEGI code (Book C cover)"
    label variable cvegeo "INEGI municipality key: state (2 digits) + municipality (3 digits)"
    keep folio mpio cvegeo
    tempfile temp_mun_c
    save `temp_mun_c'
restore
merge m:1 folio using `temp_mun_c', keep(master match) nogen

* Rural vs urban locality, from estrato in the Book C cover. We code rural = 1 for below 2,500 (estrato 4) which is INEGI's standard rural cutoff. This is a household variable with no skip pattern, the 3 rows with no folio have a missing estrato and are dropped, so the 3 households with no Book C cover stay missing, not 0
preserve
    use "$path/hh09dta_bc/c_portad.dta", clear
    keep folio estrato
    collapse (max) estrato, by(folio)
    gen rural = (estrato == 4) if !missing(estrato)
    label variable rural "Rural locality (indicator)" 
    keep folio rural
    tempfile temp_rural
    save `temp_rural'
restore
merge m:1 folio using `temp_rural', keep(master match) nogen

*******************************************************************
* Summary stats
*******************************************************************
* Restrict to children, matching MxFLS's own Book 5 (<15) scope (no obs when we do this)
*keep if age < 15
describe age gender height_cm weight_kg_gr blood_pressure_systole blood_pressure_diastole pulse hemoglobin_g_dl cholesterol cholesterol_hdl_good glycosylated_hemoglobin remittance_binary val_aid_donations_gifts moved_us us_doc child_weight_kg child_height_cm waist_cm hip_cm
summarize age gender height_cm weight_kg_gr blood_pressure_systole blood_pressure_diastole pulse hemoglobin_g_dl cholesterol cholesterol_hdl_good glycosylated_hemoglobin remittance_binary val_aid_donations_gifts moved_us us_doc child_weight_kg child_height_cm waist_cm hip_cm
* Summary stats for the new age precision / z-score variables
describe age_days age_months_exact haz waz bmiz
summarize age_days age_months_exact haz waz bmiz

* Correlation Matrix
pwcorr child_height_cm remittance_binary age, sig star(0.05)

* Testing Regressions
regress child_height_cm remittance_binary, robust
regress child_height_cm remittance_binary age, robust
regress child_weight_kg remittance_binary, robust
regress child_weight_kg remittance_binary age, robust
* Test regressions on the new HAZ/WAZ/BMIZ outcomes
regress haz remittance_binary, robust
regress haz remittance_binary age, robust
regress waz remittance_binary, robust
regress waz remittance_binary age, robust
regress bmiz remittance_binary, robust
regress bmiz remittance_binary age, robust

*******************************************************************
* Panel Identifiers & Save
*******************************************************************

* This code works. The letters being skipped 7,8 denote the wave and member type. To match the 10 format of the other waves this code does so
capture confirm string variable pid_link
if _rc tostring pid_link, replace format(%012.0f)
gen wave = 2009

* We basically grab the digits 1-6, 9-12 for individuals in waves 1 and 2. But if this is a new individual only in this wave we keep the whole thing and they get a new pid_link_canonical so they don't mistakenly try to match with a previous wave.  
gen pid_link_canonical = cond(substr(pid_link,7,1) == "C", "09_" + pid_link, substr(pid_link,1,6) + substr(pid_link,9,2) + substr(pid_link,11,2))

* Here we have unmatch 2009 with the rest of the surveys so we create a unique key for each of these individuals
replace pid_link_canonical = "09_" + folio + "_" + ls if missing(pid_link_canonical)

* Check
isid pid_link_canonical

* Same idea as 2005 but now we filling in for 2005 and 2009. 
foreach w in 2005 2002 {
    preserve
        use pid_link_canonical birth_order oldest_child birth_weight_kg delivered_home predominant_breastfed still_breastfed iron_pregnancy premature using "$clean_mxfls/mxfls_`w'_master.dta", clear
        drop if pid_link_canonical == ""
        bysort pid_link_canonical: keep if _N == 1
        foreach v in birth_order oldest_child birth_weight_kg delivered_home predominant_breastfed still_breastfed iron_pregnancy premature {
            rename `v' `v'_`w'
        }
        tempfile temp_bookiv_`w'
        save `temp_bookiv_`w''
    restore
    merge 1:1 pid_link_canonical using `temp_bookiv_`w'', keep(master match) nogen
}
foreach v in birth_order oldest_child birth_weight_kg delivered_home iron_pregnancy premature {
    replace `v' = `v'_2005 if !missing(`v'_2005)
    replace `v' = `v'_2002 if missing(`v'_2005) & !missing(`v'_2002)
}
replace predominant_breastfed = predominant_breastfed_2005 if !missing(predominant_breastfed_2005) & (still_breastfed_2005 == 0 | missing(predominant_breastfed))
replace predominant_breastfed = predominant_breastfed_2002 if !missing(predominant_breastfed_2002) & (still_breastfed_2002 == 0 | missing(predominant_breastfed))
drop *_2005 *_2002

* Household fixed effects variable. Note: That we need to grab the length, some are correct 8 number lengths and the others are not and have the lettering imbetween similar to the above. Thus we keep the correct ones and replace the ones that have the letters imbetween. Note we don't need pid_link, 11 becasue we don't need ls aka individual code
gen folio_len = length(folio)
gen folio_canonical = folio if folio_len == 8
replace folio_canonical = substr(folio,1,6) + substr(folio,9,2) if folio_len == 10

* NOTE: This folio is basically same ideas as above but for households now instead of individuals
* These cases end up with the same id bc of new individuals in wave 3 aka round c. 
gen folio_round = substr(folio,7,1) if folio_len == 10 // grabs the wave round
bysort folio_canonical: egen has_b_folio = max(folio_round == "B")
* When the folio B collides with C we have B keep its form and we give C the new folio that also contains its letters
replace folio_canonical = "09_" + folio if folio_round == "C" & has_b_folio == 1
drop folio_round has_b_folio folio_len

* CHECK stops the do file if two different 2009 households share one folio_canonical
bysort folio_canonical (folio): assert folio == folio[1]

* double checking it's unique
* isid pid_link_canonical wave

save "$clean_mxfls/mxfls_2009_master.dta", replace

*******************************************************************
* Appending Waves into Panel & Test Regressions
*******************************************************************

*** CAN CHECK THE ATTRITION ETC, BUT THE ID KEY IS GOOD, WE CHECKED.
*** Appending the two datasets on top of eachother
use "$clean_mxfls/mxfls_2009_master.dta", clear
append using "$clean_mxfls/mxfls_2005_master.dta"
append using "$clean_mxfls/mxfls_2002_master.dta"
save "$clean_mxfls/mxfls_panel_2002_2005_2009.dta", replace

* Running a couple quick test regressions
regress child_weight_kg remittance_binary, robust
regress child_weight_kg remittance_binary age, robust
regress child_height_cm remittance_binary, robust
regress child_height_cm remittance_binary age, robust

* Test regressions on the new HAZ/WAZ/BMIZ outcomes
regress haz remittance_binary, robust
regress waz remittance_binary, robust
regress bmiz remittance_binary, robust


**testing
keep if age < 15
regress hemoglobin_g_dl remittance_binary, robust
regress hemoglobin_g_dl remittance_binary age, robust



