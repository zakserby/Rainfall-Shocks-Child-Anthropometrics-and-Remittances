*******************************************************************
* Regression tables in LaTeX for the 2009 drought DID
* Project: MS Thesis
*******************************************************************

* This file exports the regression results to LaTeX the same way Summary Stats DID.do exports the summary statistics.
* It reruns the binary DID regressions of DID_Regs.do and the subanalysis do files, and reads the continuous DID results from the csv files Contdid_regs.R and the subanalysis R scripts save

* Broadly how we organize the results
	* 1) Tables/Main/regression_tables_main.tex: binary DID, SPI drought classes, continuous DID and the remittance buffer
	* 2) Tables/Appendix/regression_tables_appendix.tex: children in all three waves, subanalyses by sex and rural vs urban, event studies, continuous DID subanalyses

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
global out_sumstats  "$root/Output/Tables"
global out_regs      "$root/Output/Regression_Estimates"
global out_logs      "$root/Output/Logs"
* Create the folders this file writes to. Note: capture skips quietly if a folder already exists
foreach d in "$root/Output" "$out_sumstats" "$out_sumstats/Main" "$out_sumstats/Appendix" "$out_logs" "$out_logs/Main" {
    capture mkdir "`d'"
}

* If the master already opened the log, the capture fails quietly and the master keeps logging. If this file is run on its own, it opens the same log here and closes it at the end. Same as DID_Regs.do
capture log using "$out_logs/Main/regression_stats.log", replace text name(regstats)
local own_log = (_rc == 0)

*******************************************************************
* Some preliminaries
*******************************************************************

capture which estpost
if _rc ssc install estout
capture which distinct
if _rc ssc install distinct
capture which reghdfe
if _rc ssc install reghdfe
capture which ftools
if _rc ssc install ftools
capture which require
if _rc ssc install require

*******************************************************************
* Load DID analysis data
*******************************************************************

use "$clean_mxfls/mxfls_did_analysis.dta", clear

* Same variables the regression do files make before running
* Event study leads and lags, same as DID_Regs.do: 2002 = -2 is the lead, 2005 = -1 is the reference, 2009 = 0
gen lead_2_treat = treat2009 * (year == 2002)
gen lag_0_treat  = treat2009 * (year == 2009)
gen lead_2_remit = remittance_hh_02 * (year == 2002)
gen lag_0_remit  = remittance_hh_02 * (year == 2009)
gen lead_2_ddd   = treat2009 * remittance_hh_02 * (year == 2002)
gen lag_0_ddd    = treat2009 * remittance_hh_02 * (year == 2009)

* Children measured in all three waves, same as DID_Regs.do
foreach v in haz waz bmiz {
    bysort pid_link_canonical: egen n3_`v' = total(e_`v' == 1 & !missing(`v'))
    gen e_`v'_all3 = (e_`v' == 1 & !missing(`v') & n3_`v' == 3)
}

* Rural from the child's 2005 household, same as Subanalysis_Rural_DID_Regs.do
gen tmp05 = rural if year == 2005
bysort pid_link_canonical: egen rural_05 = max(tmp05)
drop tmp05

*******************************************************************
* Fixed Effects & Controls, same as DID_Regs.do
*******************************************************************

global FE_month "pid_link_canonical int_ym"
global did_controls "drought2005_post drought2006_post drought2007_post drought2008_post"
global controls02 "dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post rural_02_post"
global controls05 "child_male_05_post mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post"
global controls05_birth "father_height_cm_05_post father_edu_years_05_post father_age_at_birth_05_post father_employed_05_post birth_order_05_post oldest_child_05_post birth_weight_kg_05_post delivered_home_05_post predominant_breastfed_05_post iron_pregnancy_05_post premature_05_post"
* The subanalyses take out the control that does not vary within the subsample, same as the subanalysis do files
global controls02_rural "dwell_electricity_02_post dwell_floor_good_02_post dwell_roof_good_02_post dwell_ventilation_good_02_post rooms_per_person_02_post ln_pc_expend_month_02_post"
global controls05_sex "mother_height_cm_05_post mother_edu_years_05_post mother_age_at_birth_05_post mother_employed_05_post fully_vaccinated_05_post"

*******************************************************************
* Programs
*******************************************************************

* This block defines the programs used below. NOTE: Nothing runs here, just the templates. This is mainly for efficiency

* add_rows is a small helper command aka program that does this in one line, so it isn't repeated for every model, outcome and subgroup. We call it right after each reghdfe, e.g. add_rows Y Y Y N
* NOTE: We essentially are adding the number of children and 2005 municipalities, and Y or N for each set of controls
capture program drop add_rows
program define add_rows
    args droughts hh child birth
    quietly distinct pid_link_canonical if e(sample)
    local nchild = r(ndistinct)
    estadd scalar nchild = `nchild'
    estadd scalar nclust = e(N_clust)
    estadd local droughts "`droughts'"
    estadd local hh "`hh'"
    estadd local child "`child'"
    estadd local birth "`birth'"
    estadd local fe_child "Y"
    estadd local fe_month "Y"
end

* Create program that runs Models 1 to 5 of DID_Regs.do for one outcome and one sample, and stores them as m1 to m5
capture program drop run_did
program define run_did
    args v cond c02 c05
    eststo clear
    * Model 1: Basic DID
    quietly reghdfe `v' treat_post if s_`v' == 1 & `cond', absorb($FE_month) vce(cluster cvegeo05)
    add_rows N N N N
    eststo m1
    * Model 2: DID controlling for the other droughts in the window
    quietly reghdfe `v' treat_post $did_controls if s_`v' == 1 & `cond', absorb($FE_month) vce(cluster cvegeo05)
    add_rows Y N N N
    eststo m2
    * Model 3: DID with baseline 2002 household controls and 2005 child and mother controls
    quietly reghdfe `v' treat_post $did_controls `c02' `c05' if s_`v' == 1 & `cond', absorb($FE_month) vce(cluster cvegeo05)
    add_rows Y Y Y N
    eststo m3
    * Model 4: Model 3 with the father and birth controls, smaller sample
    quietly reghdfe `v' treat_post $did_controls `c02' `c05' $controls05_birth if s_`v' == 1 & `cond', absorb($FE_month) vce(cluster cvegeo05)
    add_rows Y Y Y Y
    eststo m4
    * Model 5: Triple difference, does 2002 remittance receipt buffer the drought
    quietly reghdfe `v' treat_post remittance_hh_02_post treat_post_remit02 $did_controls `c02' `c05' if s_`v' == 1 & `cond', absorb($FE_month) vce(cluster cvegeo05)
    add_rows Y Y Y N
    eststo m5
end


* Here the continuous DID results are estimated in R and saved as csv files, one row per estimate. Each row has the coefficient, its standard error, and sometimes p-values. Meaning they are not run in stata, that is why we are calling them from the CSV
* csv_cell takes one row of that csv and turns it into the three strings a cell needs:
	* 1) r(b): the coefficient with stars, e.g. -0.206**
	* 2) r(se): the standard error in parentheses, e.g. (0.089)
	* 3) r(p): the wild bootstrap p-value in brackets, e.g. [0.040]
capture program drop csv_cell
program define csv_cell, rclass
    syntax [if], se(varname) [pc(varname) wild(varname)]
	* This counts how many csv rows match the condition. There should be exactly one. If there are none that means estimate doesn't exist and it exists, doesn't guess
    quietly count `if'
    if r(N) != 1 {
        return local b ""
        return local se ""
        return local p ""
        exit
    }
    quietly summarize coef `if', meanonly
	* Coefficient
    local b = r(mean)
    quietly summarize `se' `if', meanonly
	* Standard error
    local s = r(mean)
	* This is the default two-sided p-value from a z-test, which needs to be calculated because ATT and ACRT don't come with a p-value. Note: The formula line runs for every cell, but its result is only kept for some of them
    local p = 2 * normal(-abs(`b' / `s'))
    if "`pc'" != "" {
        quietly summarize `pc' `if', meanonly
		* Here the line overwrites if there is an existing p-value and then we don't use the calculated p-value above
        local p = r(mean)
    }
    if "`wild'" != "" {
        quietly summarize `wild' `if', meanonly
        local pw = r(mean)
		* This is star rule. Stars are based on the larger of the two p-values e.g. a coefficient gets ** only if both are below .05
        local p = max(`p', `pw')
        return local p = "[" + strtrim(string(`pw', "%9.3f")) + "]"
    }
    local st ""
	* Gives the correct stars
    if `p' < 0.10 local st "\sym{*}"
    if `p' < 0.05 local st "\sym{**}"
    if `p' < 0.01 local st "\sym{***}"
	* Calling code for coefficients
    return local b = strtrim(string(`b', "%9.3f")) + "`st'"
    return local se = "(" + strtrim(string(`s', "%9.3f")) + ")"
end


* Program mostly used for counts, like the "Children" rows. It reads the number from the one row you pick, formats it with commas, and returns it as text to write into the table
capture program drop csv_num
program define csv_num, rclass
    syntax varname [if], [fmt(string)]
    if "`fmt'" == "" local fmt "%12.0fc"
    quietly count `if'
    if r(N) != 1 {
        return local n ""
        exit
    }
    quietly summarize `varlist' `if', meanonly
    return local n = strtrim(string(r(mean), "`fmt'"))
end

* Writes the start of a table that is written with file write, same LaTeX as the esttab tables
capture program drop tab_head
program define tab_head
    args ncol caption
    file write tab "\clearpage" _n "\begingroup" _n "\footnotesize" _n "\renewcommand{\arraystretch}{0.86}" _n "\setlength{\tabcolsep}{4pt}" _n "\setlength{\LTcapwidth}{\textwidth}" _n
    file write tab "\begin{longtable}{l*{`ncol'}{c}}" _n "\caption{`caption'}\\" _n "\toprule" _n
end

*******************************************************************
* LaTeX pieces used by every table
*******************************************************************

* LaTeX document setup
local preamble `" "\documentclass[12pt]{article}" "\usepackage[margin=1in]{geometry}" "\usepackage{booktabs}" "\usepackage{longtable}" "\usepackage{caption}" "\DeclareCaptionLabelSeparator{emdash}{---}" "\captionsetup{font=normalsize,labelfont=sc,textfont=sc,labelsep=emdash}" "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}" "\begin{document}" "'
local tabstart `" "\begingroup" "\footnotesize" "\renewcommand{\arraystretch}{0.86}" "\setlength{\tabcolsep}{4pt}" "\setlength{\LTcapwidth}{\textwidth}" "'
local stars "\sym{*}\(p<.10\), \sym{**}\(p<.05\), \sym{***}\(p<.01\)."

* Outcome names and age cutoffs
local cap_haz "Height-for-Age"
local cap_waz "Weight-for-Age"
local cap_bmiz "BMI-for-Age"
local head_haz "Height-for-age z-score (HAZ)"
local head_waz "Weight-for-age z-score (WAZ)"
local head_bmiz "BMI-for-age z-score (BMIZ)"
* The same names in lower case for the notes and the column headers
local low_haz "height-for-age z-score (HAZ)"
local low_waz "weight-for-age z-score (WAZ)"
local low_bmiz "BMI-for-age z-score (BMIZ)"
local age_haz "15 or younger"
local age_waz "10 or younger"
local age_bmiz "15 or younger"

* Sources cited in the notes. Used multiple times across tables so we make them locals
local cgbs "Callaway et al. (2024)"
local rural_note "Rural is a locality with fewer than 2,500 inhabitants (INEGI, 2024)."
local mnw_note "Wild cluster bootstrap p-values are in brackets, which tend to be conservative with few treated municipalities (MacKinnon et al., 2023)."
local star_note "Stars use the clustered p-values, and where a wild bootstrap p-value is shown both must pass:"

*******************************************************************
* Binary DID tables: main sample, children in all three waves, and the subanalyses
*******************************************************************

* Main idea: This blocks writes every binary DID regression table. That is 18 in all including 6 samples x 3 outcomes. 
	* The outer loop goes over the six samples: main, all3 (children in all three waves), boys, girls, rural, urban
	* The inner loop goes over the three outcomes: HAZ, WAZ, BMIZ
foreach t in main all3 boys girls rural urban {
    * Default settings. Note: The default output file is the appendix file and only main changes this
    local file "$out_sumstats/Appendix/regression_tables_appendix.tex"
    local c02 "$controls02"
    local c05 "$controls05"
	* local cond "1 == 1" puts an always-true condition in the sample-restriction spot, so the main tables use the whole estimation sample while the same run_did command still works for the subgroups
    local cond "1 == 1"
    local hh_note "2002 household controls (electricity, floor, roof, ventilation, sleeping rooms per person, log expenditure per capita, rural)"
    local child_note "2005 child and mother controls (child sex and full vaccination, mother's height, schooling, age at birth and employment)"
    local extra_note "`rural_note'"
	* This is the part of the loop that sets up each sample's settings. Before it, every table gets the main-sample defaults. Each if block then changes only what's different for that sample. i.e. 1) the sample condition, 2) the controls, 3) the title, 4) the note wording
    * The appendix tables get a shorter note that gives the reason for each check, the rest is in the main table notes. NOTE: The reason for this is because there are a bunch of tables that are very similar in terms of defining the variables but they differ only on the subset of variables e.g. boys, girls, rural, urban and so forth
    local ctrl_diff ""
    if "`t'" == "main" {
        local file "$out_sumstats/Main/regression_tables_main.tex"
        local tcap "Binary DID"
        local sample_note "The sample is children with the outcome measured in both MxFLS-2 (2005) and MxFLS-3 (2009)"
    }
    if "`t'" == "all3" {
        local cond "e_VAR_all3 == 1"
        local tcap "Binary DID, All Three Waves"
        local sample_note "The sample is children with the outcome measured in MxFLS-1 (2002), MxFLS-2 (2005) and MxFLS-3 (2009), so the event study lead and the DID come from the same children, and the regressions use 2005 and 2009"
        local reason "This robustness check keeps only children measured in all three waves (2002, 2005 and 2009), so the event study pre-trend test and the 2005--2009 DID come from the same children, and changes in sample composition, such as attrition or children first measured in 2005, cannot drive the estimates."
        local sample_short "these children"
    }
    if "`t'" == "boys" | "`t'" == "girls" {
        local c05 "$controls05_sex"
        local cond = cond("`t'" == "boys", "child_male_05 == 1", "child_male_05 == 0")
        local tcap = "Binary DID, " + cond("`t'" == "boys", "Boys", "Girls")
        local sample_note = "The sample is " + cond("`t'" == "boys", "boys", "girls") + " (2005 roster sex) with the outcome measured in both MxFLS-2 (2005) and MxFLS-3 (2009)"
        local reason "Estimating the effect separately for boys and girls tests for gender disparities in the impact of the drought, as households facing an income shock may protect the nutrition and health of sons and daughters differently."
        local sample_short = cond("`t'" == "boys", "boys", "girls") + " (2005 roster sex)"
        local ctrl_diff ", without child sex in the controls"
        local child_note "2005 mother controls and child full vaccination (mother's height, schooling, age at birth and employment), without child sex as it does not vary within the sample"
    }
    if "`t'" == "rural" | "`t'" == "urban" {
        local c02 "$controls02_rural"
        local cond = cond("`t'" == "rural", "rural_05 == 1", "rural_05 == 0")
        local tcap = "Binary DID, " + cond("`t'" == "rural", "Rural", "Urban") + " Children"
        local sample_note = "The sample is " + cond("`t'" == "rural", "rural", "urban") + " children (2005 household locality) with the outcome measured in both MxFLS-2 (2005) and MxFLS-3 (2009)"
        local reason "Estimating the effect separately for rural and urban children tests whether the drought hits harder where household incomes depend more directly on rain-fed agriculture, while urban households are mainly exposed through food prices. `rural_note'"
        local sample_short = cond("`t'" == "rural", "rural", "urban") + " children (2005 household locality)"
        local ctrl_diff ", without rural in the controls"
        local hh_note "2002 household controls (electricity, floor, roof, ventilation, sleeping rooms per person, log expenditure per capita), without rural as it barely varies within the sample"
    }

    foreach v in haz waz bmiz {
        local cond_v = subinstr("`cond'", "VAR", "`v'", .)
		* This runs Models 1–5 for this outcome
        run_did `v' "`cond_v'" "`c02'" "`c05'"
        * Building the top of the LaTeX document
        local ph `" `tabstart' "\begin{longtable}{l*{5}{c}}" "\caption{2009 Drought and `cap_`v'', `tcap'}\\" "\toprule" " & \multicolumn{5}{c}{`head_`v''}\\" "\cmidrule(lr){2-6}" "'
        local first = ("`v'" == "haz" & ("`t'" == "main" | "`t'" == "all3"))
        if `first' local ph `" `preamble' `ph' "'
        else local ph `" "\clearpage" `ph' "'
        local mode = cond(`first', "replace", "append")

        * Notes
        local pf `" "\bottomrule" "\addlinespace[6pt]\multicolumn{6}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} This table reports the effect of the 2009 drought on the WHO `low_`v'' with the binary DID. Regressions are at the child-wave level. `sample_note'," "aged `age_`v'' in both rounds. Children have a 2005 municipality and were not measured in MxFLS-3 before November 2009, as the July--October 2009 rainfall deficit had not yet fully occurred. Dropping them from both groups equalizes treated and control measurement timing. Drought equals 1 if annual CHIRPS rainfall in 2009 in the child's 2005 municipality was below its 1981--2004 20th percentile," "and post equals 1 in MxFLS-3. Remittances equal 1 if the household received cash in the last 12 months from parents, siblings, or children living outside the household in MxFLS-1 (2002), a proxy for remittances. Cols. 1--4 report the average drought" "effect. Col. 1 is the basic DID, Col. 2 adds the 2005--2008 droughts, and Col. 3 adds `hh_note' and `child_note'. Col. 4 adds the 2005 father and birth controls (father's height, schooling, age at birth and employment, birth order, oldest child," "birth weight, delivered at home, weeks predominantly breastfed, iron in pregnancy, premature), which are missing for many children, often because the father migrated. All controls are interacted with post. Col. 5 is the triple difference: it keeps the controls of Col. 3" "and adds 2002 remittances \ensuremath{\times} post and drought \ensuremath{\times} post \ensuremath{\times} 2002 remittances, the difference in the drought effect for children in 2002 remittance households. `extra_note' All regressions include child and interview year-month" "fixed effects, children in singleton cells are dropped, and standard errors are clustered by 2005 municipality, where the drought is assigned. `stars'}\\" "\end{longtable}" "\endgroup" "'
        if "`t'" != "main" local pf `" "\bottomrule" "\addlinespace[6pt]\multicolumn{6}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} `reason' The table reports the effect of the 2009 drought on the WHO `low_`v'' with the binary DID for `sample_short', aged `age_`v'' in both rounds, with a 2005 municipality and not measured in MxFLS-3 before November 2009, as the July--October 2009 rainfall deficit had not yet fully occurred. Dropping them from both groups equalizes treated and control measurement timing." "Columns, variables and controls are as in the main binary DID tables`ctrl_diff'. All regressions include child and interview year-month fixed effects, and standard errors are clustered by 2005 municipality. `stars'}\\" "\end{longtable}" "\endgroup" "'

		* This is the line that actually writes the table. It takes the five saved regressions (m1–m5) and writes them into the .tex file as five side-by-side columns. It uses the header (ph) and notes (pf) built just before, and the Y/N rows and children count that add_rows attached to each model
        esttab m1 m2 m3 m4 m5 using "`file'", `mode' booktabs ///
            cells(b(star fmt(3)) se(par fmt(3))) star(* 0.10 ** 0.05 *** 0.01) ///
            keep(treat_post treat_post_remit02 remittance_hh_02_post) order(treat_post treat_post_remit02 remittance_hh_02_post) ///
            varlabels(treat_post "Drought \ensuremath{\times} post" treat_post_remit02 "Drought \ensuremath{\times} post \ensuremath{\times} 2002 remittances" remittance_hh_02_post "2002 remittances \ensuremath{\times} post") ///
            stats(droughts hh child birth fe_child fe_month nchild, labels("2005--2008 droughts \ensuremath{\times} post" "2002 household controls \ensuremath{\times} post" "2005 child and mother controls \ensuremath{\times} post" "2005 father and birth controls \ensuremath{\times} post" "Child FE" "Interview year-month FE" "Children") fmt(0 0 0 0 0 0 %9.0fc)) ///
            collabels(none) nomtitles numbers nonotes noobs nodepvars ///
            prehead(`ph') posthead("\midrule") prefoot("\midrule") postfoot(`pf')
    }
}

*******************************************************************
* Event study tables, appendix
*******************************************************************

* Here we write the event study tables in the appendix, one per outcome. Each table has two columns: Col. 1, event study and Col. 2, triple difference event study
* 1st chunk is running the event study and DDD event study and storing it
foreach v in haz waz bmiz {
    eststo clear
    quietly reghdfe `v' lead_2_treat lag_0_treat if e_`v' == 1, absorb($FE_month) vce(cluster cvegeo05)
    quietly test lead_2_treat
    local pl = r(p)
    add_rows N N N N
    estadd scalar p_lead = `pl'
    eststo es
    quietly reghdfe `v' lead_2_ddd lag_0_ddd lead_2_treat lag_0_treat lead_2_remit lag_0_remit if e_`v' == 1 & !missing(remittance_hh_02), absorb($FE_month) vce(cluster cvegeo05)
    quietly test lead_2_ddd
    local pl = r(p)
    quietly lincom lag_0_treat + lag_0_ddd
    local lb = r(estimate)
    local ls = r(se)
    local lp = r(p)
    add_rows N N N N
    estadd scalar p_lead = `pl'
    estadd scalar lin_b = `lb'
    estadd scalar lin_se = `ls'
    estadd scalar lin_p = `lp'
    eststo ddd

	* 2nd chunk is building the header and the notes
    local ph `" "\clearpage" `tabstart' "\begin{longtable}{l*{2}{c}}" "\caption{Event Study, 2009 Drought and `cap_`v''}\\" "\toprule" " & \multicolumn{2}{c}{`head_`v''}\\" "\cmidrule(lr){2-3}" "'
    local pf `" "\bottomrule" "\addlinespace[6pt]\multicolumn{3}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} This table tests for pre-trends before the 2009 drought. Regressions are at the child-wave level and use MxFLS-1 (2002), MxFLS-2 (2005) and MxFLS-3 (2009), with 2005 as the" "reference wave. The sample is children aged `age_`v'' in each wave they are in, with a 2005 municipality, and not measured in MxFLS-3 before November 2009, as the July--October 2009 rainfall deficit had not yet fully occurred. Dropping them from both groups equalizes treated and control measurement timing. Col. 1 regresses the outcome on the drought times the 2002 wave (the lead) and times the 2009" "wave. Col. 2 adds the same terms times 2002 remittances, and the 2002 remittance terms on their own (not shown), for children with 2002 remittances observed. The P-value of the lead tests the drought lead in Col. 1 and the triple difference lead in Col. 2." "The effect for remittance households is the sum of the two 2009 drought terms. All regressions include child and interview year-month fixed effects, and standard errors are clustered by 2005 municipality. `stars'}\\" "\end{longtable}" "\endgroup" "'

	* Here we actually writing the table
    esttab es ddd using "$out_sumstats/Appendix/regression_tables_appendix.tex", append booktabs ///
        cells(b(star fmt(3)) se(par fmt(3))) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(lead_2_treat lag_0_treat lead_2_ddd lag_0_ddd) order(lead_2_treat lag_0_treat lead_2_ddd lag_0_ddd) ///
        varlabels(lead_2_treat "Drought \ensuremath{\times} 2002 wave (lead)" lag_0_treat "Drought \ensuremath{\times} 2009 wave" lead_2_ddd "Drought \ensuremath{\times} 2002 remittances \ensuremath{\times} 2002 wave (lead)" lag_0_ddd "Drought \ensuremath{\times} 2002 remittances \ensuremath{\times} 2009 wave") ///
        stats(fe_child fe_month p_lead lin_b lin_se lin_p N nchild nclust, labels("Child FE" "Interview year-month FE" "P-value, lead = 0" "\midrule 2009 drought effect, remittance households" "\quad SE" "\quad P-value" "\midrule Observations (child-waves)" "Children" "2005 municipalities (clusters)") fmt(0 0 3 3 3 3 %9.0fc %9.0fc %9.0fc)) ///
        collabels(none) mtitles("Event study" "Triple difference") numbers nonotes noobs nodepvars ///
        prehead(`ph') posthead("\midrule") prefoot("\midrule") postfoot(`pf')
}

*******************************************************************
* Load the continuous DID results from R
*******************************************************************

* Each csv is imported as text and appended, then the numbers are turned back into numbers. Note: NA in the csv files becomes missing
tempfile csvs
local first 1
foreach f in "Main/contdid_results" "Main/contdid_spi_results" "Appendix/subanalysis_gender_contdid_results" "Appendix/subanalysis_gender_contdid_spi_results" "Appendix/subanalysis_rural_contdid_results" "Appendix/subanalysis_rural_contdid_spi_results" {
    import delimited using "$out_regs/`f'.csv", clear varnames(1) stringcols(_all) bindquote(strict)
    gen file = "`f'"
    if `first' == 0 append using `csvs'
    save `csvs', replace
    local first 0
}
* This turns each number column from text back into numbers
foreach v in coef se_contdid_bootstrap se_clustered_municipality n_children se_clustered p_clustered p_wild_bootstrap children_in_group municipalities_in_group untreated_municipalities {
    capture destring `v', replace force
}

* The stacked csv files label their rows in different columns, so we make two columns that work for every row. That is spec (the label) and group (all or the subgroup)
* Then each table cell is found with one condition, e.g. outcome == "haz" & group == "boys" & spec == "..."

* The SPI files put the label in specification and the continuous DID files put it in estimate, so spec takes whichever one the row has
gen spec = specification
replace spec = estimate if spec == ""
replace term = "" if missing(term)
* group is "all" for the main results. Subgroup labels end in ", boys", ", girls", ", rural" or ", urban", in spec or in sample
gen group = "all"
foreach g in boys girls rural urban {
    * If the label contains ", boys" (or girls, rural, urban), that row belongs to that subgroup
    replace group = "`g'" if strpos(spec, ", `g'") | strpos(sample, ", `g'")
}

* For assumptions the continuous DID and SPI notes use the same sentence 
local pt_note "The average ATT needs parallel trends, while the average ACRT, the average slope of the dose response, and comparisons between doses need strong parallel trends (Callaway et al., 2024)."

*******************************************************************
* SPI drought classes table, main
*******************************************************************

* NOTE: This writes Table 7, the SPI drought-class table, into the Main tables file. The SPI regressions were run in R, so esttab can't be used. Instead the table is written by hand, one LaTeX line at a time with file write and using numbers pulled from the csv rows with csv_cell and csv_num
* NOTE: Columns 1, 3 and 5 are dry (z <= -1) vs near normal, columns 2, 4 and 6 split dry into moderately dry and severely dry or worse, from section 8 of Contdid_regs.R

* Creating the header
capture file close tab
file open tab using "$out_sumstats/Main/regression_tables_main.tex", write append
tab_head 6 "2009 Drought and Child Anthropometrics, SPI Drought Classes"
file write tab " & \multicolumn{2}{c}{HAZ} & \multicolumn{2}{c}{WAZ} & \multicolumn{2}{c}{BMIZ}\\" _n "\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}" _n " & (1) & (2) & (3) & (4) & (5) & (6)\\" _n "\midrule" _n

* Coefficient rows: each outcome has the dry column and the groups column
foreach r in dry mod sev {
	* This loops over the three row types: dry, moderately dry, and severely dry or worse. And ensure math prints the math 
    if "`r'" == "dry" local lab "Dry, \ensuremath{z\le -1}"
    if "`r'" == "mod" local lab "Moderately dry, \ensuremath{-1.5<z\le -1}"
    if "`r'" == "sev" local lab "Severely dry or worse, \ensuremath{z\le -1.5}"
	* Three rows the coefficient with the label, the SEs and the p-value from the wild bootstrap
    local row1 "`lab'"
    local row2 ""
    local row3 ""
    foreach v in haz waz bmiz {
        foreach c in dry groups {
			* Picking the row 
            local cnd `"group == "all" & outcome == "`v'""'
			* adds the dry-vs-near-normal regression
            if "`c'" == "dry" local cnd `"`cnd' & strpos(spec, "Average ATT: dry") & term == "dry" & "`r'" == "dry""'
			* For a groups column, this adds the split regression, the term mod_dry or sev_dry
            if "`c'" == "groups" local cnd `"`cnd' & strpos(spec, "Groups: moderately") & term == "`r'_dry""'
			* Calls the row and gives the se, p-value, and wild boostrap p value
            csv_cell if `cnd', se(se_clustered) pc(p_clustered) wild(p_wild_bootstrap)
            local row1 "`row1' & `r(b)'"
            local row2 "`row2' & `r(se)'"
            local row3 "`row3' & `r(p)'"
        }
    }
    file write tab "`row1'\\" _n "`row2'\\" _n "`row3'\\" _n
}

* Pre-trend p-values (wild bootstrap), number of children and municipalities
* The midrule just creates that medium thick line
file write tab "\midrule" _n 
foreach r in mod sev {
	* This just labels the rows 
    local row = "P-val: pre-trend 2002--2005, " + cond("`r'" == "mod", "moderately dry", "severely dry or worse")
    foreach v in haz waz bmiz {
		* Grabs wild bootstrap p-values and formats it 
        csv_num p_wild_bootstrap if group == "all" & outcome == "`v'" & strpos(spec, "Pre-trend") & term == "`r'_dry", fmt(%9.3f)
		* The first & skips the first column as the event study is only for the split cases not the dry column
        local row "`row' & & `r(n)'"
    }
	* This writes the row
    file write tab "`row'\\" _n
}
* CSV holds the count for each of these three rows, each with a label, we just call upon each 
foreach n in n_children municipalities_in_group untreated_municipalities {
    if "`n'" == "n_children" local row "Children"
    if "`n'" == "municipalities_in_group" local row "Dry 2005 municipalities"
    if "`n'" == "untreated_municipalities" local row "Near-normal 2005 municipalities"
    foreach v in haz waz bmiz {
        csv_num `n' if group == "all" & outcome == "`v'" & strpos(spec, "Average ATT: dry") & term == "dry"
        local a "`r(n)'"
        csv_num `n' if group == "all" & outcome == "`v'" & strpos(spec, "Groups: moderately") & term == "mod_dry"
        local b "`r(n)'"
        if "`n'" == "municipalities_in_group" local b ""
		* Adds the two cells for each outcome
        local row "`row' & `a' & `b'"
    }
	* Writes the row
    file write tab "`row'\\" _n
}
* Notes
file write tab "\bottomrule" _n
file write tab "\addlinespace[6pt]\multicolumn{7}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} This table compares children in dry 2005 municipalities with children in near-normal ones. Regressions are at the child level, on the 2005 to 2009 change in the WHO z-score," _n
file write tab "for the binary DID estimation samples without children in wet municipalities (z of 1 or above). z is 2009 annual CHIRPS rainfall in the child's 2005 municipality in standard deviations from its 1981--2004 mean, grouped with the SPI classes of" _n
file write tab "McKee et al. (1993): near normal is z between -1 and 1, moderately dry between -1.5 and -1, and severely dry or worse -1.5 or below. Cols. 1, 3 and 5 compare dry (z of -1 or below) with near normal, and Cols. 2, 4 and 6 split dry into moderately dry and" _n
file write tab "severely dry or worse. Each class is only compared with near normal, which needs parallel trends, and not with the other class, which would need strong parallel trends (Callaway et al., 2024). Standard errors clustered by 2005 municipality are in parentheses." _n
file write tab "`mnw_note' `star_note' `stars'}\\" _n
file write tab "\end{longtable}" _n "\endgroup" _n

*******************************************************************
* Continuous DID and remittance buffer table, main
*******************************************************************

* This writes the table, "Continuous DID and Remittances as a Moderator". Like Table 7, it's written by hand with file write, because the results come from R csv files.
tab_head 3 "Continuous DID and Remittances as a Moderator"
file write tab " & HAZ & WAZ & BMIZ\\" _n " & (1) & (2) & (3)\\" _n "\midrule" _n

* Writing the header
local rows ""
* Adding the rows to panel A, this has 4 rows and is seperated by "|"
local rows `"`rows' "A|\textit{Panel A: Continuous DID, any 2009 rainfall deficit}" "c|Average ATT|spi|Main average ATT|D_any" "c|Average ACRT|res|Average ACRT|-" "n|Children|spi|Main average ATT|D_any""'
* Adding the rows to panel B and is seperated by "|"
local rows `"`rows' "A|\addlinespace[6pt]\textit{Panel B: Average ATT by 2002 remittances}" "c|No 2002 remittances|res|2002 remittances = 0|-" "c|2002 remittances|res|2002 remittances = 1|-" "c|Difference, remittances minus none|res|Average ATT difference|-" "n|Children|res|Average ATT difference|-""'
* Adding the rows to panel C and is seperated by "|"
local rows `"`rows' "A|\addlinespace[6pt]\textit{Panel C: Dry vs near normal by 2002 remittances}" "c|Dry, \ensuremath{z\le -1}|spi|Remittance buffer|dry" "c|Dry \ensuremath{\times} 2002 remittances|spi|Remittance buffer|dry:remittance_hh_02" "c|2002 remittances|spi|Remittance buffer|remittance_hh_02" "n|Children|spi|Remittance buffer|dry""'
foreach r of local rows {
	* This cuts the text at each "|" and numbers the pieces. The "|" signs are kept as pieces too
	* NOTE: (defined as) 1 = row type, 3 = label, 5 = csv file type, 7 = specification, 9 = term
    tokenize "`r'", parse("|")
    if "`1'" == "A" {
        file write tab "`3' & & & \\" _n
        continue
    }
	* The contdid csv has no term, so - in the row list becomes a blank term
    local tm = cond("`9'" == "-", "", "`9'")
	* Again this the coefficient, SE, wild P-val
    local row1 "`3'"
    local row2 ""
    local row3 ""
	* For each outcome, this builds the condition that finds one csv row
    foreach v in haz waz bmiz {
		* Describes the row we want
        local cnd `"group == "all" & outcome == "`v'" & strpos(spec, "`7'") & term == "`tm'""'
        if "`1'" == "n" {
			* So it looks up the children count for this outcome and adds it to the table row, this is for counts
            csv_num n_children if `cnd'
            local row1 "`row1' & `r(n)'"
        }
		* Finds the matching row and returns the three things listed below, this is for coefficients
        else if "`5'" == "spi" {
            csv_cell if `cnd', se(se_clustered) pc(p_clustered) wild(p_wild_bootstrap)
            local row1 "`row1' & `r(b)'"
            local row2 "`row2' & `r(se)'"
            local row3 "`row3' & `r(p)'"
        }
		* This works the same way, but uses our municipality bootstrap SE, and there's no wild p-value, this is for coefficents
        else {
            csv_cell if `cnd', se(se_clustered_municipality)
            local row1 "`row1' & `r(b)'"
            local row2 "`row2' & `r(se)'"
        }
    }
	* Writing it
    file write tab "`row1'\\" _n
    if "`1'" == "c" file write tab "`row2'\\" _n
    if "`1'" == "c" & "`5'" == "spi" file write tab "`row3'\\" _n
}
* Line
file write tab "\bottomrule" _n
* Notes
file write tab "\addlinespace[6pt]\multicolumn{4}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} This table reports the continuous DID of `cgbs' and tests whether 2002 remittances moderated the effect of the 2009 drought. Regressions are at the child level, on the 2005 to 2009" _n
file write tab "change in the WHO z-score, for the binary DID estimation samples. The dose is the 2009 rainfall deficit, max(0,~-z), where z is 2009 annual CHIRPS rainfall in the child's 2005 municipality in standard deviations from its 1981--2004 mean. Panel A compares" _n
file write tab "children with any 2009 deficit (z below 0) with children without one, Panel B gives the average ATT separately for children in 2002 remittance and non-remittance households and their difference, and Panel C interacts dry (z of -1 or below) vs near normal" _n
file write tab "with 2002 remittances. Remittances are cash received in the last 12 months from parents, siblings, or children living outside the household in MxFLS-1 (2002), a proxy for remittances. `pt_note' The estimator does not allow controls or fixed effects." _n
file write tab "Standard errors in parentheses are clustered by 2005 municipality, where the drought is assigned. Because the package only bootstraps at the child level, we coded a multiplier bootstrap of its influence functions that clusters by 2005 municipality, and the other rows use municipality-clustered regression standard errors. `mnw_note' `star_note' `stars'}\\" _n
file write tab "\end{longtable}" _n "\endgroup" _n "\end{document}" _n
file close tab

*******************************************************************
* Continuous DID by subgroup tables, appendix
*******************************************************************

* This writes the continuous DID subgroup tables in the appendix, one per outcome. It's the same approach as Tables 7 and 8, meaning it's hand-written with file write, with the numbers pulled from the csv data with csv_cell and csv_num
* It's slightly different as the columns are the subgroups (boys, girls, rural, urban) and the rows are the estimates, this is to safe space and make less tables

* Opening the appendix file 
file open tab using "$out_sumstats/Appendix/regression_tables_appendix.tex", write append
foreach v in haz waz bmiz {
	* This writes the table's opening
    tab_head 4 "Continuous DID by Subgroup, `cap_`v''"
	* Writes column headings one at a time
    file write tab " & \multicolumn{4}{c}{Change in the `low_`v''}\\" _n "\cmidrule(lr){2-5}" _n " & Boys & Girls & Rural & Urban\\" _n " & (1) & (2) & (3) & (4)\\" _n "\midrule" _n
    foreach r in any dry sev remit {
		* Here we kinda looping through stuff and then giving them the full labels with the mathematical symbol as well
        if "`r'" == "any" local lab "Any 2009 deficit"
        if "`r'" == "dry" local lab "Dry vs near normal, \ensuremath{z\le -1}"
        if "`r'" == "sev" local lab "Severely dry or worse vs near normal, \ensuremath{z\le -1.5}"
        if "`r'" == "remit" local lab "Dry \ensuremath{\times} 2002 remittances"
		* Again this the coefficient, SE, wild P-val
        local row1 "`lab'"
        local row2 ""
        local row3 ""
        foreach g in boys girls rural urban {
			* Finding csv row
            local cnd `"group == "`g'" & outcome == "`v'""'
			* Each line adds the regression label and term for this estimate, so the condition picks exactly one csv row
            if "`r'" == "any" local cnd `"`cnd' & strpos(spec, "Main average ATT") & term == "D_any""'
            if "`r'" == "dry" local cnd `"`cnd' & strpos(spec, "Average ATT: dry") & term == "dry""'
            if "`r'" == "sev" local cnd `"`cnd' & strpos(spec, "Groups: moderately") & term == "sev_dry""'
            if "`r'" == "remit" local cnd `"`cnd' & strpos(spec, "Remittance buffer") & term == "dry:remittance_hh_02""'
			* This fetches the coefficient with stars, the clustered SE, and the wild p-value
            csv_cell if `cnd', se(se_clustered) pc(p_clustered) wild(p_wild_bootstrap)
			* Saving these three lines copy those numbers onto the end of row1, row2 and row3
            local row1 "`row1' & `r(b)'"
            local row2 "`row2' & `r(se)'"
            local row3 "`row3' & `r(p)'"
        }
		* Writes the three finished lines, then moves to the next estimate
        file write tab "`row1'\\" _n "`row2'\\" _n "`row3'\\" _n
    }
	* Here we seperate out the coefficients with the sample sizes. The gist of this code is to grab the counts
    foreach n in n_children untreated_municipalities {
        local row = cond("`n'" == "n_children", "\midrule Children", "Municipalities with no 2009 deficit")
        foreach g in boys girls rural urban {
            csv_num `n' if group == "`g'" & outcome == "`v'" & strpos(spec, "Main average ATT") & term == "D_any"
            local row "`row' & `r(n)'"
        }
        file write tab "`row'\\" _n
    }
	* Line at the bottom 
    file write tab "\bottomrule" _n
	* Notes
    file write tab "\addlinespace[6pt]\multicolumn{5}{p{0.97\textwidth}}{\scriptsize \textit{Notes:} This table repeats the continuous DID and SPI results for the WHO `low_`v'' by subgroup. Cols. 1--2 test for gender disparities in the impact of the drought (2005 roster sex)," _n
    file write tab "and Cols. 3--4 test whether it hits rural children harder, as rural household incomes depend more directly on rain-fed agriculture. `rural_note' Rows, samples (aged `age_`v'' in both rounds) and definitions are as in the main continuous DID and SPI" _n
    file write tab "tables, and subgroups with few dry or untreated municipalities give noisy estimates. Standard errors clustered by 2005 municipality are in parentheses. `mnw_note' `star_note' `stars'}\\" _n
    file write tab "\end{longtable}" _n "\endgroup" _n
}
* Writing it and then end
file write tab "\end{document}" _n
file close tab


if `own_log' log close regstats






