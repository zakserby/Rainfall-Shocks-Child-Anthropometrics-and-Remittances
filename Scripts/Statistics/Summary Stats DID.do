*******************************************************************
* Summary statistics for the 2009 drought DID
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

* Create the folders this project writes to (capture skips quietly if a folder already exists)
foreach d in "$root/Data" "$root/Data/Cleaned" "$clean_mxfls" "$clean_who" "$clean_rain" "$root/Output" "$out_sumstats" "$out_regs" "$out_logs" {
    capture mkdir "`d'"
}
* This includes writing to Main and Appendix subfolders, same as 00_Master.do
foreach d in "$out_sumstats" "$out_regs" "$out_logs" {
    foreach s in Main Appendix {
        capture mkdir "`d'/`s'"
    }
}

*******************************************************************
* Some preliminaries
*******************************************************************

capture which estpost
if _rc ssc install estout
capture which distinct
if _rc ssc install distinct

*******************************************************************
* Load DID analysis data
*******************************************************************

* Loading in data
use "$clean_mxfls/mxfls_did_analysis.dta", clear

*******************************************************************
* DID variables
*******************************************************************

* Keep remittances only on the 2005 row, blank on the others
gen remit05_tmp = remittance_hh if year == 2005
* Copy that 2005 value to all of the child's rows
bysort pid_link_canonical: egen remittance_hh_05 = max(remit05_tmp)
drop remit05_tmp

gen did_sample = age < 16 & !missing(treat2009) & treat_early == 0

*******************************************************************
* Globals
*******************************************************************

* Controls. We sort sections by household, child and mother, father and birth. Primarily bc father and birth have fewer observations for children and reduce the sample size a lot 
* Outcomes
global outcomes "haz waz bmiz"
* Treatment and the other droughts between the waves
global treatment "treat2009 drought2005 drought2006 drought2007 drought2008"
* Remittances
global remit "remittance_hh_05"
* Household controls
global controls_hh "dwell_electricity dwell_floor_good dwell_roof_good dwell_ventilation_good rooms_per_person ln_pc_expend_month rural"
* Child and mother controls
global controls_child "age child_male fully_vaccinated mother_height_cm mother_edu_years mother_age_at_birth mother_employed"
* Father and Book IV birth controls. Note this is filled for fewer children
global controls_birth "father_height_cm father_edu_years father_age_at_birth father_employed birth_order oldest_child birth_weight_kg delivered_home predominant_breastfed iron_pregnancy premature"

* All vars
global baseline_vars "$outcomes $remit $controls_hh $controls_child $controls_birth"
* Table variables as they enter the regressions, so 2002 remittances and household controls, and 2005 child and parent controls
global remit02 "remittance_hh_02"
global controls_hh02 "dwell_electricity_02 dwell_floor_good_02 dwell_roof_good_02 dwell_ventilation_good_02 rooms_per_person_02 ln_pc_expend_month_02 rural_02"
* Continuous DID treatment
global contdid_treat "rain_deficit2009 any_deficit2009"
global table_vars "$treatment $contdid_treat $remit02 $controls_hh02 $controls_child $controls_birth"


*******************************************************************
* Summary statistics table, 2005 at baseline
*******************************************************************

* Creating general summary table for our primary regression samples. Split into two panels. 
local cells "mean(fmt(2)) sd(fmt(2)) min(fmt(2)) p50(fmt(2)) max(fmt(2)) count(fmt(%9.0fc))"
local groups `"treat2009 "\addlinespace[6pt]\textit{Binary DID treatment (CHIRPS)}" rain_deficit2009 "\addlinespace[6pt]\textit{Continuous DID treatment (CHIRPS)}" remittance_hh_02 "\addlinespace[6pt]\textit{Remittances (MxFLS-1, 2002)}" dwell_electricity_02 "\addlinespace[6pt]\textit{Household (MxFLS-1, 2002)}" age "\addlinespace[6pt]\textit{Child and mother (MxFLS-2, 2005)}" father_height_cm "\addlinespace[6pt]\textit{Father and birth (MxFLS-1 and MxFLS-2)}""'

* We compute the numbers in the notes underneath the table here, so nothing is hardcoded. Same regressions as the sample size table further down including Models 1 and 3 of DID_Regs.do
capture which reghdfe
if _rc ssc install reghdfe
* reghdfe needs ftools and require, so on a new computer they must be installed before the first reghdfe below
capture which ftools
if _rc ssc install ftools
capture which require
if _rc ssc install require
foreach o in haz waz bmiz {
    quietly count if reg_`o' == 1 & year == 2005
	* n_ = children in the estimation sample
    local n = r(N)
    quietly reghdfe `o' treat_post if s_`o' == 1, absorb(pid_link_canonical int_ym) vce(cluster cvegeo05)
    quietly distinct pid_link_canonical if e(sample)
	* single_ = children Model 1 drops as singletons
    local single_`o' = `n' - r(ndistinct)
    quietly reghdfe `o' treat_post drought2005_post drought2006_post drought2007_post drought2008_post dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post rural_02_post child_male_05_post mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post if s_`o' == 1, absorb(pid_link_canonical int_ym) vce(cluster cvegeo05)
    quietly distinct pid_link_canonical if e(sample)
	* Here we add some formatting for these 4 lines below
	* m3_ = children in Model 3
    local m3_`o' : display %9.0fc r(ndistinct)
    local m3_`o' = strtrim("`m3_`o''")
    quietly count if reg_`o' == 1 & year == 2005 & rain_z2009 >= 1 & !missing(rain_z2009)
	* wet_ = children with z of 1 or above that the SPI analysis drops
    local wet_`o' = r(N)
    local n_`o' : display %9.0fc `n'
    local n_`o' = strtrim("`n_`o''")
}

* Panel A: WAZ sample
preserve
    keep if reg_waz == 1 & year == 2005
    quietly estpost summarize waz $table_vars, detail
    esttab using "$out_sumstats/Main/summary_stats_did_table.tex", replace booktabs label ///
        cells("`cells'") collabels("Mean" "Std. dev." "Min." "Median" "Max." "Obs.") noobs nonumber nomtitle ///
        refcat(waz "\textit{Outcomes (MxFLS-2, 2005)}" `groups', nolabel) ///
        prehead("\documentclass[12pt]{article}" "\usepackage[margin=1in]{geometry}" "\usepackage{booktabs}" "\usepackage{longtable}" "\usepackage{caption}" ///
            "\DeclareCaptionLabelSeparator{emdash}{---}" "\captionsetup{font=normalsize,labelfont=sc,textfont=sc,labelsep=emdash}" "\begin{document}" ///
            "\begingroup" "\footnotesize" "\renewcommand{\arraystretch}{0.86}" "\setlength{\tabcolsep}{4pt}" "\setlength{\LTcapwidth}{\textwidth}" ///
            "\begin{longtable}{l*{6}{c}}" "\caption{Summary Statistics, 2009 Drought DID Samples}\\" ///
            "\toprule" " & \multicolumn{6}{c}{Panel A: WAZ, age \ensuremath{\le}10 in 2005 and 2009}\\" "\cmidrule(lr){2-7}") ///
        posthead("\midrule") prefoot("") postfoot("\bottomrule" "\newpage")
restore

* Panel B: HAZ and BMIZ sample, appended to the same table
preserve
    keep if reg_haz == 1 & year == 2005
    quietly estpost summarize haz bmiz $table_vars, detail
    esttab using "$out_sumstats/Main/summary_stats_did_table.tex", append booktabs label ///
        cells("`cells'") collabels("Mean" "Std. dev." "Min." "Median" "Max." "Obs.") noobs nonumber nomtitle ///
        refcat(haz "\textit{Outcomes (MxFLS-2, 2005)}" `groups', nolabel) ///
        prehead("\toprule" " & \multicolumn{6}{c}{Panel B: HAZ and BMIZ, age \ensuremath{\le}15 in 2005 and 2009}\\" "\cmidrule(lr){2-7}") ///
        posthead("\midrule") prefoot("") ///
        postfoot("\bottomrule" ///
            "\addlinespace[6pt]\multicolumn{7}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} Sample is children with the outcome measured in MxFLS-2 (2005) and MxFLS-3 (2009), aged 10 or younger (Panel A) or 15 or younger (Panel B)" ///
            "in both rounds, with a 2005 municipality and not measured in MxFLS-3 before November 2009. Each row is a child, and Obs. counts children. Statistics are at the 2005 baseline, except" ///
            "remittances and household variables (2002). Panel B reports BMIZ for the HAZ sample, as the BMIZ sample is similar and slightly smaller, at `n_bmiz' children." ///
            "The continuous DID uses these samples, and binary DID Model 1 excludes `single_waz' (Panel A) and `single_haz' (Panel B) children in singleton" ///
            "interview-month cells. Model 3, the main binary DID, adds 2005--2008 droughts and 2002 household and 2005 child and mother controls and includes" ///
            "`m3_waz' (WAZ), `m3_haz' (HAZ), and `m3_bmiz' (BMIZ) children. The SPI analysis compares z of -1 or below with z between -1 and 1, excluding" ///
            "`wet_waz' (Panel A) and `wet_haz' (Panel B) children with z of 1 or above. Drought is annual CHIRPS rainfall in the 2005 municipality below its 1981--2004 20th percentile." ///
            "Rainfall deficit is max(0,~-z), where z is 2009 rainfall in standard deviations from the 1981--2004 municipal mean. Remittances, a proxy, are cash received in the last" ///
            "12 months from parents, siblings, or children outside the household. Expenditure is in 2002 pesos per person per month, a proxy for" ///
            "household income. Rural is a locality under 2,500 inhabitants (INEGI, 2024).}\\" ///
            "\end{longtable}" "\endgroup")
restore


*******************************************************************
* Balance table of treated vs control
*******************************************************************


* reghdfe needed for absorbing fixed effects; install once if not already present, same as in DID_Regs
capture which reghdfe
if _rc ssc install reghdfe
capture which ftools
if _rc ssc install ftools
capture which require
if _rc ssc install require


* This block defines a new Stata command called cnt_row. NOTE: Nothing runs here, just the template. This program writes one row of the table by 
	* 1) Counting how many control children, households or municipalities meet a condition
	* 1) Counting how many treated children, households or municipalities meet a condition
	* 3) Writes one line of LaTeX into the .tex file that label control count & treated count.
* Has about 60 rows, and this is done because more efficient then writing 60 times. Only the label and conditions change
capture program drop cnt_row
program define cnt_row
    gettoken label 0 : 0
    syntax , Sample(string) Control(string) Treated(string) [Distinct(varname) Last]
    foreach g in control treated {
		* Counting Children
		* NOTE: Basically box left blank no distinct(...) in the line then the program sees an empty box and counts children, one per row of data. 
        if "`distinct'" == "" {
            quietly count if (`sample') & (``g'')
            local n_`g' = r(N)
        }
		* Counting municipalities or households
		* NOTE: This distinct(...) is filled in with cvegeo05 or folio_canonical  here and specified later down
        else {
            capture quietly distinct `distinct' if (`sample') & (``g'')
            local n_`g' = cond(_rc, 0, r(ndistinct))
        }
		* Just formatting
        local c_`g' : display %12.0fc `n_`g''
        local c_`g' = strtrim("`c_`g''")
    }
    local end "\\*"
	* LaTeX stuff
    if "`last'" != "" local end "\\"
    file write cnt `"`label' & `c_control' & `c_treated'`end'"' _n
	* Displaying in stata
    di as text `"`label': "' as result "`c_control' control, `c_treated' treated"
end


preserve
* Run DID Models 1, 3 and 5 of DID_Regs.do, we flag with e(sample), meaning this creates a row with 1's and 0's and the 1's are the observations used in the regression. Then later on in the script we will count these ones to get the sample of treated or control. 
* Purpose of this block is to find out exactly which children each regression uses
local cnt_did_controls "drought2005_post drought2006_post drought2007_post drought2008_post"
local cnt_controls02 "dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post rural_02_post"
local cnt_controls05 "child_male_05_post mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post"
foreach o in haz waz bmiz {
    quietly reghdfe `o' treat_post if s_`o' == 1, absorb(pid_link_canonical int_ym) vce(cluster cvegeo05)
    gen m1_`o' = e(sample)
    quietly reghdfe `o' treat_post `cnt_did_controls' `cnt_controls02' `cnt_controls05' if s_`o' == 1, absorb(pid_link_canonical int_ym) vce(cluster cvegeo05)
    gen m3_`o' = e(sample)
    quietly reghdfe `o' treat_post remittance_hh_02_post treat_post_remit02 `cnt_did_controls' `cnt_controls02' `cnt_controls05' if s_`o' == 1, absorb(pid_link_canonical int_ym) vce(cluster cvegeo05)
    gen m5_`o' = e(sample)
}

* Every child in the regression samples has exactly two rows used, 2005 and 2009. The table counts children, not rows, so each child must be counted once, thus we only keep one year
keep if year == 2005
* SPI classes from z = 2009 rainfall in the 2005 municipality, the same cutoffs as add_spi() in functions.R
gen spi_near = (rain_z2009 > -1) if rain_z2009 < 1
gen spi_mod = (rain_z2009 <= -1 & rain_z2009 > -1.5) if rain_z2009 < 1
gen spi_sev = (rain_z2009 <= -1.5) if rain_z2009 < 1
gen spi_sevonly = (rain_z2009 <= -1.5 & rain_z2009 > -2) if rain_z2009 < 1
gen spi_ext = (rain_z2009 <= -2) if rain_z2009 < 1

* Text for the group headings
local age_haz "15 or younger"
local age_waz "10 or younger"
local age_bmiz "15 or younger"

capture file close cnt
file open cnt using "$out_sumstats/Main/summary_stats_did_table.tex", write append 

* LaTex
file write cnt "\clearpage" _n "\begingroup" _n "\footnotesize" _n "\renewcommand{\arraystretch}{0.86}" _n "\setlength{\tabcolsep}{4pt}" _n "\setlength{\LTcapwidth}{\textwidth}" _n "\setlength{\LTleft}{0pt}" _n "\setlength{\LTright}{0pt plus 1fil}" _n
file write cnt "\begin{longtable}{@{\extracolsep{\fill}}lcc}" _n "\caption{Sample Sizes, Treated and Control Children}\\" _n 


*******************************************************************
* Panel A: binary DID
*******************************************************************
* LaTeX
file write cnt "\toprule" _n " & \multicolumn{2}{c}{Panel A: Binary DID, 2009 drought}\\*" _n "\cmidrule(lr){2-3}" _n " & Control & Treated\\*" _n "\midrule" _n
foreach o in haz waz bmiz {
    local O = upper("`o'")
    local space ""
    if "`o'" != "haz" local space "\addlinespace[6pt]"
    file write cnt "`space'\textit{`O' sample, children aged `age_`o'' in 2005 and 2009} & & \\*" _n
    local t "control(treat2009 == 0) treated(treat2009 == 1)"
	* This is calling upon the count if in our program and the t sets control vs treated counts
    cnt_row "Children, estimation sample", sample(reg_`o' == 1) `t'
    cnt_row "Children, Model 1 (basic DID)", sample(m1_`o' == 1) `t'
    cnt_row "Children, Model 3 (DID with controls)", sample(m3_`o' == 1) `t'
    cnt_row "Children, Model 5 (triple difference)", sample(m5_`o' == 1) `t'
	* Only exists for triple difference, and also adds up to total number of children in Model 5 with the children not in remittance households
    cnt_row "\quad In 2002 remittance households", sample(m5_`o' == 1 & remittance_hh_02 == 1) `t'
    cnt_row "Households, estimation sample", sample(reg_`o' == 1) `t' distinct(folio_canonical)
    cnt_row "2005 municipalities, estimation sample", sample(reg_`o' == 1) `t' distinct(cvegeo05) last
}
file write cnt "\bottomrule" _n

*******************************************************************
* Panel B: continuous DID, any 2009 rainfall deficit vs none
*******************************************************************

file write cnt "\addlinespace[12pt]" _n "\toprule" _n " & \multicolumn{2}{c}{Panel B: Continuous DID, any 2009 deficit}\\*" _n "\cmidrule(lr){2-3}" _n " & Control & Treated\\*" _n "\midrule" _n
foreach o in haz waz bmiz {
    local O = upper("`o'")
    local space ""
    if "`o'" != "haz" local space "\addlinespace[6pt]"
    file write cnt "`space'\textit{`O' sample, children aged `age_`o'' in 2005 and 2009} & & \\*" _n
    local t "control(any_deficit2009 == 0) treated(any_deficit2009 == 1)"
    cnt_row "Children, continuous DID sample", sample(reg_`o' == 1) `t'
	* Does not add up to total number of households thus we have both variables
    cnt_row "\quad In 2002 remittance households", sample(reg_`o' == 1 & remittance_hh_02 == 1) `t'
    cnt_row "\quad In 2002 non-remittance households", sample(reg_`o' == 1 & remittance_hh_02 == 0) `t'
    cnt_row "Households", sample(reg_`o' == 1) `t' distinct(folio_canonical)
    cnt_row "2005 municipalities", sample(reg_`o' == 1) `t' distinct(cvegeo05) last
}
file write cnt "\bottomrule" _n

*******************************************************************
* Panel C: continuous DID, SPI drought classes vs near normal
*******************************************************************

* Panel C starts on a new page like Panel B of the summary statistics table, so its header is not left at the bottom of a page
file write cnt "\newpage" _n "\toprule" _n " & \multicolumn{2}{c}{Panel C: SPI classes vs near normal}\\*" _n "\cmidrule(lr){2-3}" _n " & Control & Treated\\*" _n "\midrule" _n
* Each class with its z range, in the order of the SPI regressions in Contdid_regs.R
local cls_mod "Moderately dry, \ensuremath{-1.5<z\le -1}"
local cls_sev "Severely dry or worse, \ensuremath{z\le -1.5}"
local cls_sevonly "Severely dry, \ensuremath{-2<z\le -1.5}"
local cls_ext "Extremely dry, \ensuremath{z\le -2}"
foreach c in mod sev sevonly ext {
    local space ""
    if "`c'" != "mod" local space "\addlinespace[6pt]"
    file write cnt `"`space'\textit{`cls_`c''} & & \\*"' _n
    local t "control(spi_near == 1) treated(spi_`c' == 1)"
    cnt_row "Children, HAZ sample", sample(reg_haz == 1) `t'
    cnt_row "Children, WAZ sample", sample(reg_waz == 1) `t'
    cnt_row "Children, BMIZ sample", sample(reg_bmiz == 1) `t'
    cnt_row "2005 municipalities, HAZ sample", sample(reg_haz == 1) `t' distinct(cvegeo05)
    cnt_row "2005 municipalities, WAZ sample", sample(reg_waz == 1) `t' distinct(cvegeo05)
    cnt_row "2005 municipalities, BMIZ sample", sample(reg_bmiz == 1) `t' distinct(cvegeo05) last
}

* Notes at the bottom
file write cnt "\bottomrule" _n
file write cnt "\addlinespace[6pt]\multicolumn{3}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} Counts are children, households, or 2005 municipalities in each regression sample, at the 2005 baseline." _n
file write cnt "The estimation samples are children with the outcome measured in both MxFLS-2 (2005) and MxFLS-3 (2009), aged 15 or younger (HAZ and BMIZ) or 10 or younger (WAZ) in both rounds." _n
file write cnt "Children have a 2005 municipality and were not measured in MxFLS-3 before November 2009. In Panel A, treated children's 2005 municipality had a drought in 2009, with annual CHIRPS" _n 
file write cnt "rainfall below its 1981--2004 20th percentile. Model 1 drops children alone in their interview-month cell, and the continuous DID, with no interview-month fixed effects, uses the full sample. Model 3 adds the 2005--2008 droughts and 2002 household" _n 
file write cnt "and 2005 child and mother controls, and keeps children with all of them observed. Model 5 adds 2002 remittances times post and the drought times post times 2002 remittances, and keeps children" _n
file write cnt "with 2002 remittances observed, so Model 5 children not in remittance households are the Model 5 count minus the remittance row. In Panel B, treated children's 2005 municipality had 2009" _n
file write cnt "rainfall below its 1981--2004 mean (z below 0), where z is 2009 rainfall in standard deviations from that mean. These are the continuous DID samples, and the two remittance rows are the groups" _n
file write cnt "of the remittance split. They do not add up to the total, as children with missing 2002 remittances are in neither row. Panel C compares each SPI class with near normal (z between -1 and 1)" _n
file write cnt "in the continuous DID. As McKee, Doesken and Kleist (1993) note, drought intensity is arbitrarily defined for values of the SPI with these categories. Children in wet municipalities (z of 1 or above)" _n
file write cnt "are not in the SPI analysis. Remittances are cash received in the last 12 months" _n
file write cnt "from parents, siblings, or children living outside the household, from MxFLS-1 (2002).}\\" _n
file write cnt "\end{longtable}" _n "\endgroup" _n

file close cnt
restore



*******************************************************************
* Attrition
*******************************************************************

* Table 3: Sample attrition across the waves, as counts of children in each wave. The DID only uses 2005 and 2009 and the event studies use all three waves. The continuous DID uses the same samples, so this covers it too
* Each row starts from the children who meet the regression sample rules in the first wave and could stay in the sample, not from all children measured. Then the later columns count who is still in the sample

preserve
* Age cutoffs of the samples, same as in Build_DID_Analysis_Data.do
local cut_haz 16
local cut_waz 11
local cut_bmiz 16

* The loop copies each child's age in 2002, 2005 and 2009 onto all of that child's rows. That way, each row can compare ages across waves, for example age2009-age2005. A missing age2009 marks the children whose 2009 age has to be filled in some other way: from their birth date and the 2009 interview date, or failing that, from the median gap.
foreach y in 2002 2005 2009 {
    gen tmp = age if year == `y'
    bysort pid_link_canonical: egen age`y' = max(tmp)
    drop tmp
}

* Median age gap between the waves, for children under 16 measured in both
* We use the median bc a few gaps are weird, like −4 or 44, probably from data entry errors
gen gap0509 = age2009 - age2005 if year == 2005 & age < 16
quietly summarize gap0509, detail
local gap0509 = r(p50)
local n0509 = r(N)
* Share of these children within a year of the median, for the notes section
quietly count if inrange(gap0509, `gap0509' - 1, `gap0509' + 1)
local sh0509 : display %3.0f 100 * r(N) / `n0509'
local sh0509 = strtrim("`sh0509'")
gen gap0209 = age2009 - age2002 if year == 2002 & age < 16
quietly summarize gap0209, detail
local gap0209 = r(p50)
local n0209 = r(N)
quietly count if inrange(gap0209, `gap0209' - 1, `gap0209' + 1)
local sh0209 : display %3.0f 100 * r(N) / `n0209'
local sh0209 = strtrim("`sh0209'")

* Copy each child's 2009 interview date to all of their rows
gen tmp = int_date if year == 2009
bysort pid_link_canonical: egen int2009 = max(tmp)
drop tmp

* This code below has essentially three steps
* 1) Measured: if the child was measured in 2009, use that real age.
* 2) Computed: if they're in MxFLS-3 but weren't measured, compute their age from their birth date and their household's 2009 interview date. This only happens if the birth date passed the cleaning check. This code is similar to the cleaning script
* 3) Estimated: if they're not in MxFLS-3 (or the birth date failed the check), take their last measured age plus the median gap: 4 years from 2005, 8 from 2002.
gen age_next = age2009
replace age_next = floor((int2009 - mdy(bmonth, 15, byear)) / 365.25) if missing(age_next) & !missing(age_days)
replace age_next = age2005 + `gap0509' if missing(age_next)
replace age_next = age2002 + `gap0209' if missing(age_next)

capture file close cnt
file open cnt using "$out_sumstats/Main/summary_stats_did_table.tex", write append

* LaTeX
file write cnt "\clearpage" _n "\begingroup" _n "\footnotesize" _n "\renewcommand{\arraystretch}{0.86}" _n "\setlength{\tabcolsep}{4pt}" _n "\setlength{\LTcapwidth}{\textwidth}" _n "\setlength{\LTleft}{0pt}" _n "\setlength{\LTright}{0pt plus 1fil}" _n
file write cnt "\begin{longtable}{@{\extracolsep{\fill}}lccc}" _n "\caption{Sample Attrition Across Waves}\\" _n
file write cnt "\toprule" _n " & \multicolumn{3}{c}{Number of children (N)}\\*" _n "\cmidrule(lr){2-4}" _n " & MxFLS-1 (2002) & MxFLS-2 (2005) & MxFLS-3 (2009)\\*" _n "\midrule" _n

foreach o in haz waz bmiz {
    local O = upper("`o'")
    * This codes children as 1 if children meet the DID sample rules for the starting baseline regression specification
    gen did = (year == 2005 & age < `cut_`o'' & !missing(`o') & !missing(treat2009) & treat_early == 0 & age_next < `cut_`o'')
    * If children meet the event study rules
    gen es = (year == 2002 & age < `cut_`o'' & !missing(`o') & !missing(cvegeo) & treat_early == 0 & age_next < `cut_`o'')
    bysort pid_link_canonical: egen es05 = max(e_`o' == 1 & !missing(`o') & year == 2005)
    bysort pid_link_canonical: egen es_n3 = total(e_`o' == 1 & !missing(`o'))
    * Counting children
    quietly count if did == 1
    local d05 = r(N)
    quietly count if did == 1 & reg_`o' == 1
    local d09 = r(N)
    quietly count if es == 1
    local e02 = r(N)
    quietly count if es == 1 & es05 == 1
    local e05 = r(N)
    quietly count if es == 1 & es_n3 == 3
    local e09 = r(N)
    * Just formatting
    foreach x in d05 d09 e02 e05 e09 {
        local `x' : display %12.0fc ``x''
        local `x' = strtrim("``x''")
    }
    * LaTeX stuff
    local space ""
    if "`o'" != "haz" local space "\addlinespace[6pt]"
    file write cnt "`space'\textit{`O' sample, children under `cut_`o''} & & & \\*" _n
    file write cnt "DID (2005 and 2009) & & `d05' & `d09'\\*" _n
    file write cnt "Event studies (all three waves) & `e02' & `e05' & `e09'\\" _n
    * Displaying in stata
    di as text "`O': DID 2005 `d05', 2009 `d09' | event studies 2002 `e02', 2005 `e05', 2009 `e09'"
    drop did es es05 es_n3
}

* Notes at the bottom, like Tables 1 and 2
file write cnt "\bottomrule" _n
file write cnt "\addlinespace[6pt]\multicolumn{4}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} Counts of children in each wave. The DID only uses MxFLS-2 (2005) and MxFLS-3 (2009), and the event studies" _n
file write cnt "use all three waves. The continuous DID uses the same samples. Each row starts from the children who meet the regression sample rules in its first wave and could stay in the sample, not" _n
file write cnt "from all children measured. They are under the age cutoff in every wave, 15 or younger (HAZ and BMIZ) or 10 or younger (WAZ), have a municipality, and were not measured in MxFLS-3 before" _n 
file write cnt "November 2009. Later columns count the children still in the sample, so the drop across columns is attrition. This is not the same as the children the regressions drop as singletons, with" _n
file write cnt "only one wave in the sample, as these also include children who aged out or are only in a later wave. The DID 2009 count is the estimation sample, and the event study 2009 count is the" _n
file write cnt "children in all three waves. For children not measured in a later wave, the age then comes from their birth date and their household's MxFLS-3 measurement date, as MxFLS-3 interviews ran" _n 
file write cnt "over several years. Children not in MxFLS-3, or whose birth date does not agree with their reported age, get their last measured age plus the median age gap between the waves, `gap0509' years" _n
file write cnt "from 2005 to 2009 and `gap0209' years from 2002 to 2009. This is an estimate, but `sh0509' and `sh0209' percent of these gaps are within a year of the median among children measured in both" _n
file write cnt "waves, and a guess off by a year only matters for children near the age cutoff, who could be counted as attrition when they aged out or the other way around. So those who would have aged out anyway are not counted as attrition.}\\" _n
file write cnt "\end{longtable}" _n "\endgroup" _n

* The document ends here
file write cnt "\end{document}" _n
file close cnt
restore










