/*************************************************************************************************************/
/********  	WHO 2007 Reference                      	                                  	  ****/
/********  	Department of Nutrition for Health and Development                                 	  ****/
/********  	World Health Organization                                                          	  ****/
/********  	Last modified on 08/008/2007     -     For STATA versions 7.0 and above             	  ****/
/*************************************************************************************************************/

/*************************************************************************************************************/
/********  	A macro / program for calculating the z-scores and prevalences for a nutritional survey  *****/
/*************************************************************************************************************/

*! version 1	  08aug2007
program define who2007
	version 7.0
	discard
	set type float
	args reflib datalib datalab sex age ageunit weight height oedema sw 
	if ("`age'" == "" | "`ageunit'" == "" | "`weight'" == "" | "`height'" == "" | "`oedema'" == "" | "`sw'" == "") {
		di "Error - You must specify 10 arguments: "
		di "	reflib datalib datalab sex age ageunit weight height oedema sw"
	exit
	}
	/*	Phase I =  generation of z-scores = _zwfa, _zhfa, _zbfa */
	di as txt _n "Please wait, programme is running............."
	di as txt _n ".............................................."

	/*	Data preparation */
	
	/* quietly begins*/
	qui {
	tempvar tsort tsex toedema 
	tempvar ageint agelow agehigh agediff
	tempvar zwfa asd3p asd23p asd3n asd23n 
	tempvar zhfa zbfa bsd3p bsd23p bsd3n bsd23n 
		 
	gen `tsort'=_n
	inspect `sex'
	if  r(N_unique)==0 {
		gen `tsex'=1 if `sex'=="m" | `sex'=="M" 
		replace `tsex'=2 if `sex'=="f" | `sex'=="F"
		replace `tsex'=. if `sex'== " " 
	}
	if  r(N_unique)~=0 {
		gen `tsex'=1 if `sex'==1
		replace `tsex'=2 if `sex'==2
		replace `tsex'=. if `sex'==.
	}
	qui ta `oedema'
	if  r(r)==0 {
		gen `toedema'=1  if `oedema'==.
	}
	if  r(r)~=0 {
		gen `toedema'= 1 
		replace `toedema'= 2 if `oedema'=="y" | `oedema'=="Y"
	}
	
	/**	Age calculations /interpolations **/
	capture gen double _agemons=`age' 
	if _rc==0 {
		replace _agemons=`age'/30.4375 if `ageunit'=="days"
		replace _agemons=`age'*12 if `ageunit'=="years"
		lab var _agemons "Calculated age in months for deriving z scores"	
	}
	if _rc~=0 {
		drop _agemons
		gen double _agemons=`age'
		replace _agemons=`age'/30.4375 if `ageunit'=="days"
		replace _agemons=`age'*12 if `ageunit'=="years"
		lab var _agemons "Calculated age in months for deriving z scores"	
	}
	
	tempvar agetemp 
	gen double `agetemp'=`age'
	replace `agetemp'=`age'/30.4375 if `ageunit'=="days"	
	replace `agetemp'=`age'*12 if `ageunit'=="years"
	
	macro def outage1 "if _agemons < 61 | _agemons >=121"
	macro def outage2 "if _agemons < 61 | _agemons >=229"
	
	macro def inage1 "_agemons >= 61 & _agemons <121"
	macro def inage2 "_agemons >= 61 & _agemons <229"
	
	capture gen double _cbmi= `weight'*10000/(`height'*`height')
	if _rc==0 {
		lab var _cbmi "Calculated bmi=weight / squared(height)"
	}
	if _rc~=0 {
		drop _cbmi
		gen double _cbmi= `weight'*10000/(`height'*`height')
		lab var _cbmi "Calculated bmi=weight / squared(height)"
	}		
	
	* ========================== Interpolation of l,m,s values ============================
	
	gen double `ageint'=_agemons
	gen double `agelow'=int(_agemons)
	gen double `agehigh'=`agelow'+1 
	gen double `agediff'=_agemons-`agelow'
	
	* ===============wfalow LMS calculations============
	replace `ageint'=`agelow'
	sort `tsex' `ageint'
	local string "xxx\wfawho2007.dta"
	local i=`reflib'
	global wfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$wfafile"

	foreach Y in  l m s {
		gen double `Y'1 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 
	
	* ===============wfahigh LMS calculations============
	replace `ageint'=`agehigh'
	sort `tsex' `ageint'
	local string "xxx\wfawho2007.dta"
	local i=`reflib'
	global wfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$wfafile"

	foreach Y in  l m s {
		gen double `Y'2 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 

	foreach Y in  l m s {
		gen double `Y' =`Y'1*(1+`agelow') - (`Y'2*`agelow') + (`Y'2-`Y'1)*_agemons
		*gen double `Y' =`Y'1 + `agediff'*(`Y'2-`Y'1)
	}
	gen double `zwfa'=(((`weight'/m)^l)-1)/(s*l)
	gen double `asd3p'=m*((1+l*s*3)^(1/l))
	gen double `asd23p'=`asd3p'- m*((1+l*s*2)^(1/l))
	replace `zwfa'= 3+((`weight'-`asd3p')/`asd23p') if (`zwfa'>3 & `zwfa'~=.) 
	gen double `asd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `asd23n'= m*((1+l*s*(-2))^(1/l))-`asd3n'
	replace `zwfa'=-3-((`asd3n'-`weight')/`asd23n') if (`zwfa'<-3 & `zwfa'~=.) 
	
	gen double _zwfa= `zwfa'
	replace _zwfa =. $outage1
	drop l l1 l2 m m1 m2 s s1 s2
	
	* ==========================  ZHFA  =======================================

	* ===============hfalow LMS calculations============
	replace `ageint'=`agelow'
	sort `tsex' `ageint'
	local string "xxx\hfawho2007.dta"
	local i=`reflib'
	global hfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$hfafile"

	foreach Y in  l m s {
		gen double `Y'1 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 
	
	* ===============hfahigh LMS calculations============
	replace `ageint'=`agehigh'
	sort `tsex' `ageint'
	local string "xxx\hfawho2007.dta"
	local i=`reflib'
	global hfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$hfafile"

	foreach Y in  l m s {
		gen double `Y'2 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 

	foreach Y in  l m s {
		gen double `Y' =`Y'1*(1+`agelow') - (`Y'2*`agelow') + (`Y'2-`Y'1)*_agemons
	}

	gen double `zhfa'=(((`height'/m)^l)-1)/(s*l)
	
	gen double _zhfa= `zhfa'
	replace _zhfa =. $outage2
	drop l l1 l2 m m1 m2 s s1 s2
	
	* =============================== ZBFA ======================================

	* ===============bfalow LMS calculations============
	replace `ageint'=`agelow'
	sort `tsex' `ageint'
	local string "xxx\bfawho2007.dta"
	local i=`reflib'
	global bfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$bfafile"

	foreach Y in  l m s {
		gen double `Y'1 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 
	
	* ===============bfahigh LMS calculations============
	replace `ageint'=`agehigh'
	sort `tsex' `ageint'
	local string "xxx\bfawho2007.dta"
	local i=`reflib'
	global bfafile: subinstr local string "xxx" "`i'"
	merge `tsex' `ageint' using "$bfafile"

	foreach Y in  l m s {
		gen double `Y'2 =`Y'
	}
	keep if _merge~=2
	drop l m s _merge 

	foreach Y in  l m s {
		gen double `Y' =`Y'1*(1+`agelow') - (`Y'2*`agelow') + (`Y'2-`Y'1)*_agemons
	}
	gen double `zbfa'=(((_cbmi/m)^l)-1)/(s*l)
	gen double `bsd3p'=m*((1+l*s*3)^(1/l))
	gen double `bsd23p'=`bsd3p'- m*((1+l*s*2)^(1/l))
	replace `zbfa'= 3+((_cbmi-`bsd3p')/`bsd23p') if (`zbfa'>3 & `zbfa'~=.) 
	gen double `bsd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `bsd23n'= m*((1+l*s*(-2))^(1/l))-`bsd3n'
	replace `zbfa'=(-3)-((`bsd3n'-_cbmi)/`bsd23n') if (`zbfa'<-3 & `zbfa'~=.) 
		
	gen double _zbfa= `zbfa'
	replace _zbfa =. $outage2
	drop l l1 l2 m m1 m2 s s1 s2
	
	set type float
	*===================================================================================
		
	foreach Y in  _zhfa _zwfa _zbfa {
		replace `Y' =round(`Y', 0.01)
	}

	* ====	Set weight-based z-scores to missing for Oedema cases =====
	
	foreach Y in _zwfa _zbfa {
		replace `Y'=. if `toedema'==2
	}
	
	gen _fwfa=0 if _zwfa~=.
	gen _fhfa=0 if _zhfa~=.
	gen _fbfa=0 if _zbfa~=.
	
	replace _fhfa=1 if (_zhfa< -6 | _zhfa >6)
	replace _fwfa=1 if (_zwfa< -6 | _zwfa >5)
	replace _fbfa=1 if (_zbfa< -5 | _zbfa >5)
	
	foreach Y in wfa bfa hfa {
		replace _f`Y'=. if _z`Y'==.
	}
	
	lab var _zhfa "Height-for-age z-score"
	lab var _zwfa "Weight-for-age z-score"
	lab var _zbfa "BMI-for-age z-score"
	
	lab var _fbfa "=1 if (_zbfa < -5 | _zfa >5)"
	lab var _fhfa "=1 if (_zhfa < -6 | _zhfa >6)"
	lab var _fwfa "=1 if (_zwfa < -6 | _zwfa >5)"
	
	/*	 Clean-up after Phase I = sort data as originally provided	*/
	
	sort `tsort'
	
	} /* quietly ends*/
		
	di "Note: z-scores are flagged according to the following rules:"
	di " "
	di "				_zhfa = . if (_zhfa < -6 or _zhfa >6)"
	di "				_zwfa = . if (_zwfa < -6 or _zwfa >5)"
	di "				_zbfa = . if (_zbfa < -5 or _zbfa >5)"
	di "				"
	
	drop __*
	local string "xxx\yy_z.dta"
	local i=`datalib'
	local j=`datalab'
	global outf: subinstr local string "xxx" "`i'" 
	global outfile: subinstr global outf "yy" "`j'"
	save "$outfile", replace
	di "Note 1:	Original data plus z-scores are written to"
	di "		$outfile"
	di "				"
	di "				"
	local string "xxx\yy_z.xls"
	local i=`datalib'
	local j=`datalab'
	global outs: subinstr local string "xxx" "`i'" 
	global outsheet: subinstr global outs "yy" "`j'"
	outsheet using "$outsheet", replace
	di "Note 2: 	Original data plus z-scores are written to"
	di "		$outsheet"
	
	di "				"
	di as txt _n "Please wait, programme is calculating prevalences............."
	di as txt _n ".............................................."
	
	
	/*	Phase II=	Generation of prevalences */
	
	/*	Check the sampling weights before generating prevalences */
	
	qui summ `sw'
	if r(min)<0	{
		di "Error - Negative sampling weights encountered, prevalence tables are not produced!"
	exit
	}
	
	if r(min)>=0 {
			/* quietly begins */
	qui {
	
	tempvar agetemp hcw1 noage agegrp
	
	gen `hcw1'=int(_agemons/12)
	gen `noage'=1 if `hcw1'==.
	gen `agegrp'=`hcw1'
	recode `agegrp' 5=1 6=2 7=3 8=4 9=5 10=6 11=7 12=8 13=9 14=10 15=11 16=12 17=13 18=14 19=15 20/max=.
	replace `agegrp'=. if  _agemons<61
	replace `agegrp'=. if  _agemons>=229
	replace `agegrp'=0 if `noage'==1
	lab def agegrp  0 "no age" 1 "5" 2 "6" 3 "7" 4 "8" 5 "9" 6 "10" 7 "11" 8 "12" 9 "13" 10 "14" 11 "15" 12 "16" 13 "17" 14 "18" 15 "19" , modify 
	lab val `agegrp' agegrp
	
	/*	Here I needed to declare the temporary variables for sex, oedema and ovflag again */
	
	tempvar tsex toedema ovflag
	inspect `sex'
	if  r(N_unique)==0 {
		gen `tsex'=1 if `sex'=="m" | `sex'=="M" 
		replace `tsex'=2 if `sex'=="f" | `sex'=="F"
		replace `tsex'=. if `sex'== " " 
	}
	if  r(N_unique)~=0 {
		gen `tsex'=1 if `sex'==1
		replace `tsex'=2 if `sex'==2
		replace `tsex'=. if `sex'==.
	}

	qui tab `oedema'
	if  r(r)==0 {
		gen `toedema'=1  if `oedema'==.
	}
	if  r(r)~=0 {
		gen `toedema'= 1 
		replace `toedema'= 2 if `oedema'=="y" | `oedema'=="Y"
	}
	tempvar ovflag _fhfa _fwfa _fbfa
	
	gen `_fhfa'=0 if _zhfa~=.
	gen `_fwfa'=0 if _zwfa~=.
	gen `_fbfa'=0 if _zbfa~=.
	
	replace `_fhfa'=1 if (_zhfa< -6 | _zhfa >6)
	replace `_fwfa'=1 if (_zwfa< -6 | _zwfa >5)
	replace `_fbfa'=1 if (_zbfa< -5 | _zbfa >5)
	
	gen `ovflag'=((`_fhfa'==1) | (`_fwfa'==1) | (`_fbfa'==1))
	
	foreach X in _zwfa _zhfa _zbfa  {
		tempvar `X'p1 `X'p2 `X'p3 `X'm2 `X'm3
	}
		
	/*	For each indicator, declare binary variables to calculate prevalences*/
		*	p1= 1 if zscore above 1SD,  0 otherwise
		*	p2= 1 if zscore above 2SD,  0 otherwise
		*	p3= 1 if zscore above 3SD,  0 otherwise
		*	m2= 1 if zscore below -2SD, 0 otherwise
		*	m3= 1 if zscore below -3SD, 0 otherwise
	
	foreach Y in _zhfa {
		gen ``Y'p1'= `Y' > 1 & `Y' ~=. 
		gen ``Y'p2'= `Y' > 2 & `Y' ~=. 
		gen ``Y'p3'= `Y' > 3 & `Y' ~=. 
		gen ``Y'm2'=  (`Y' < -2 & `Y' ~=.) 
		gen ``Y'm3'=  (`Y' < -3 & `Y' ~=.) 
	}
	
	/*	Special adjustment for the weight-based indicators with oedema for m2 and m3*/
	
	foreach Y in _zwfa _zbfa {
		gen ``Y'p1'= `Y' > 1 & `Y' ~=. 
		gen ``Y'p2'= `Y' > 2 & `Y' ~=. 
		gen ``Y'p3'= `Y' > 3 & `Y' ~=. 
		gen ``Y'm2'=  (`Y' < -2 & `Y' ~=.) 
		replace ``Y'm2'= 1 if `toedema'==2 
		gen ``Y'm3'=  (`Y' < -3 & `Y' ~=.) 
		replace ``Y'm3'= 1 if `toedema'==2 
	}
	
	
	set type double	
	
		/*	declare flagged values missing by indicator*/
		replace _zwfa=. if `_fwfa'==1
		replace _zhfa=. if `_fhfa'==1
		replace _zbfa=. if `_fbfa'==1
		
		/*	Sexes combined */
		
		/*	Start with the zscores for Weight-for-age	*/		
		
		foreach X in _zwfa {
		
			/*	Calculates Means and SDs - for total agegroups */
		
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=.				/* temp=1 is a tag for all children with non-missing _zwfa */
			gen `swt'=`X'*`sw' if `temp'==1 		/* swt=_zwfa * sampling weight for all children if their temp=1 */
			egen `y'=sum(`swt') if `temp'==1		/* y=add-up swt for all children if their temp=1 */
			egen `x'=sum(`sw')  if `temp'==1		/* x=add-up sampling weight of all children if their temp=1 */
			gen `yx'=(`y'/`x')				/* yx=weighted mean of _zwfa= y/x	*/
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1  	/* swtsd= sum of squared (_zwfa - mean) times sampling weight */
			egen `ysd'=sum(`swtsd') if `temp'==1		/* ysd=add-up swtsd of all children if their temp=1 */
			gen `yxsd'=(`ysd'/(`x'-1))^0.5			/* yxsd=standard deviation=ysd divided by x-1*/
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)	/* put the mean _zwfa in matrix mean0*/
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'					/* put the SD _zwfa in matrix sd0*/
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}		
			/*	Calculates Means and SDs - for disaggregated agegroups (same logic as above)*/
			
			forvalues Z = 1/6 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw')  if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}		
			}
			
			/*	Collate row matrices of total and disaggregate ages
				into matrix mean and matrix sd of weight-for-age zscore	*/
			
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 
			}
			
			/*	Calculation of prevalences starts here */
			
			foreach Y in m3 m2 p1 p2 p3 {
				
				/*	Calculates prevalences for total agegroups */
				
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=. & $inage1		 /* temp=1 is a tag for all children with non-missing _zwfa */
				recode `temp' .=1 if `toedema'==2 & $inage1  	/* also temp=1 if children have oedema & in age range*/
				
				/*	in the "recode" above, is a special adjustment for oedema cases which applies to _zwfa _zbfa*/
					 
						
				/*	Take for example Y=m3, i.e. the prevalence of _zwfa below -3SD	*/
				
				gen `swt'=``X'`Y''*`sw' if `temp'==1 		/* swt=_zwfam3 times sampling weight for all children if their temp=1*/
				egen `y'=sum(`swt') if `temp'==1		/* y=add-up swt for all children if their temp=1*/
				egen `x'=sum(`sw') if `temp'==1		/* x=add-up sampling weight of all children if their temp=1 */	
				gen `yx'=(`y'/`x')				/* yx=required prevalence of _zwfa below -3SD	*/	
				summ `x'					/* x=Weighted N = sum of sampling weights of all children if their temp=1*/
				if r(N)~=0 {
					matrix A0=r(mean)			/* put the weighted N in matrix A0*/
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)			/* put the proportion yx in matrix B0*/
				}
				if r(N)==0 {
					matrix B0= -69
				}
				
				/*	Calculates prevalences - for disaggregated agegroups (same logic as above)*/
				
				 forvalues Z = 1/6  {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & $inage1
					recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z' & $inage1
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				/* 	For total (Z=0) and disaggregate ages (Z= 1 to 6) 
					matrix A holds the weighted N 
					matrix B holds the proportion or prevalence
					matrix X holds lower 95% CI
					matrix Y holds upper 95% CI	
				*/					
				
				forvalues Z = 0/6  {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
				
				/*	Collates  row matrices of total and disaggregate ages */
				
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 
				}
				
				/*	Collates column matrices of prevalence (B) and lower 95% (X) and higher 95% (Y)*/
				matrix `Y'=B, X, Y
			}
			
			/*	Collating together the full matrix for _zwfa
				A = Weighted N (which is the same whether the prevalence is m2 or m3)
				m3 = prevalence m3 with 95%C.I.
				m2 = prevalence m2 with 95%C.I.
				mean = mean zscore weight-for-age (_zwfa)
				sd =   SD zscore weight-for-age (_zwfa)
			*/				
			matrix `X'= A, m3, m2, p1, p2, p3, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI "  "Mean" "SD" 
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
			matrix rownames `X'="Total" "5" "6" "7" "8" "9" "10"
		}
		
		/**	Height_for_age	**/
		foreach X in _zhfa {
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=. 
			gen `swt'=`X'*`sw' if `temp'==1 
			egen `y'=sum(`swt') if `temp'==1
			egen `x'=sum(`sw')  if `temp'==1
			gen `yx'=(`y'/`x')
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			egen `ysd'=sum(`swtsd') if `temp'==1
			gen `yxsd'=(`ysd'/(`x'-1))^0.5
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}
			forvalues Z = 1/15 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}
			}
		
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
			}
		
			foreach Y in m3 m2 p1 p2 p3  {
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=. & $inage2
				gen `swt'=``X'`Y''*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				summ `x'
				if r(N)~=0 {
					matrix A0=r(mean)
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)
				}
				if r(N)==0 {
					matrix B0= -69
				}
				forvalues Z = 1/15 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & $inage2
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				forvalues Z = 0/15 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
		
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
				}
				matrix `Y'=B, X, Y
			}
			matrix `X'= A, m3, m2, p1, p2, p3, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
			matrix rownames `X'="Total" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19"
		}
		
		foreach X in _zbfa {
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=. 
			gen `swt'=`X'*`sw' if `temp'==1 
			egen `y'=sum(`swt') if `temp'==1
			egen `x'=sum(`sw') if `temp'==1
			gen `yx'=(`y'/`x')
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			egen `ysd'=sum(`swtsd') if `temp'==1
			gen `yxsd'=(`ysd'/(`x'-1))^0.5
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}
			forvalues Z = 1/15 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z' 
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}
			}
		
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
			}
			foreach Y in m3 m2 p1 p2 p3 {
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=. 
				recode `temp' .=1 if `toedema'==2 & $inage2
				gen `swt'=``X'`Y''*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				summ `x'
				if r(N)~=0 {
					matrix A0=r(mean)
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)
				}
				if r(N)==0 {
					matrix B0= -69
				}
				forvalues Z = 1/15 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z'
					recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				forvalues Z = 0/15 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
		
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
				}
				matrix `Y'=B, X, Y
			}
			matrix `X'= A, m3, m2, p1,p2, p3, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
			matrix rownames `X'="Total" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19"
		}
		
		/* Sex=1  Males;	 Sex=2 Females */
		foreach S in 1 2  {
			foreach X in _zwfa {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1 
				egen `x'=sum(`sw') if `temp'==1 
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1 
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				forvalues Z = 1/6 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			 		egen `ysd'=sum(`swtsd') if `temp'==1 
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
		
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 
				}
				foreach Y in m3 m2 p1 p2 p3 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S' & $inage1
					recode `temp' .=1 if `toedema'==2 & `agegrp'~=. & `tsex' ==`S' & $inage1
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
					forvalues Z = 1/6 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S' & $inage1
						recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z' & `tsex' ==`S' & $inage1
						gen `swt'=``X'`Y''*`sw' if `temp'==1 
						egen `y'=sum(`swt') if `temp'==1 
						egen `x'=sum(`sw') if `temp'==1 
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					 forvalues Z = 0/6 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, p1,p2, p3, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
				matrix rownames `X'`S'="Total" "5" "6" "7" "8" "9" "10" 
			}
			foreach X in _zhfa {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				forvalues Z = 1/15 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
					egen `ysd'=sum(`swtsd') if `temp'==1
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
			
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
				}
			
				foreach Y in m3 m2 p1 p2 p3 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S' & $inage2
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x') 
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
					forvalues Z = 1/15 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S' & $inage2
						gen `swt'=``X'`Y''*`sw' if `temp'==1
						egen `y'=sum(`swt') if `temp'==1
						egen `x'=sum(`sw') if `temp'==1
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					forvalues Z = 0/15 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, p1, p2, p3, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
				matrix rownames `X'`S'="Total" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19"
			}
			foreach X in _zbfa {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1 
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				forvalues Z = 1/15 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
					egen `ysd'=sum(`swtsd') if `temp'==1
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
			
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
				}
				foreach Y in m3 m2 p1 p2 p3 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S' & $inage2 
					recode `temp' .=1 if `toedema'==2 & `agegrp'~=. & `tsex' ==`S' & $inage2
					gen `swt'=``X'`Y''*`sw' if `temp'==1
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
	
					forvalues Z = 1/15 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S' & $inage2
						recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z' & `tsex' ==`S' & $inage2
						gen `swt'=``X'`Y''*`sw' if `temp'==1
						egen `y'=sum(`swt') if `temp'==1
						egen `x'=sum(`sw') if `temp'==1
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					forvalues Z = 0/15 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6 \ `Q'7 \ `Q'8 \ `Q'9 \ `Q'10 \ `Q'11 \ `Q'12 \ `Q'13 \ `Q'14 \ `Q'15
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, p1,p2, p3, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
				matrix rownames `X'`S'="Total" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19"
			}
			
	} /* quietly ends*/
	
	} 	/* end of loop "if r(min)~=0" */

	/**		On-screen display	**/
			
	di as txt _n "Set 1:		Sexes combined"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwfa, noheader 
	di as txt _n "Age groups	Height-for-age"
	matrix list _zhfa, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbfa, noheader 
	
	di as txt _n "Set 2:		Males"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwfa1, noheader 
	di as txt _n "Age groups	Height-for-age"
	matrix list _zhfa1, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbfa1, noheader 
	
	di as txt _n "Set 3:		Females"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwfa2, noheader 
	di as txt _n "Age groups	Height-for-age"
	matrix list _zhfa2, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbfa2, noheader  
	
	/**	Prevalence tables to the excel sheet begin here	**/
	
	/* quietly begins */
	qui {
	#delimit ;
	gen __t1=-1000; gen __t2=-2000; gen __t3=-9000; gen __t4=-10000; gen __t5=-3000; gen __t6=-9000; 
	gen __t7=-10000; gen __t8=-6000; gen __t9=-9000; gen __t10=-10000; gen __t11=-7000; gen __t12=-9000;
	gen __t13=-10000; gen __t14=-8000; gen __t15=-9000; gen __t16=-10000; gen __t17=-4000; gen __t18=-5000;
	#delimit cr
	
	/*
	forval i = 1/7 {
		gen __u`i'=__t`i'
	}
	gen __u8=-4000
	gen __u9=-5000
	
	mkmat __u1-__u9 in 1 , matrix(U)*/
	
	mkmat __t1-__t18 in 1, matrix(T)
	
	/*matrix emp=J(7,9,-6666)
	matrix uemp=J(1,9,-6666)
	forval i = 1/2 {
		matrix _zwfa`i'=_zwfa`i', emp
	}
	matrix _zwfa=_zwfa, emp
	matrix U=U, uemp
	
	*set trace on*/
	
	matrix first=J(141,1,0)
	forval i = 1/141 {
		matrix first[`i',1]=`i'
	}
	matrix row1=J(1,18,-6666)
	matrix row1[1,1]=-11	/*	Set1:*/
	matrix row1[1,2]=-22	/*	Sexes*/
	matrix row1[1,3]=-33	/*	combined*/
	
	matrix rowM=J(1,18,-6666)
	matrix rowM[1,1]=-1155	/*	Set2:*/
	matrix rowM[1,2]=-2255	/*	Males*/
	
	matrix rowF=J(1,18,-6666)
	matrix rowF[1,1]=-1166	/*	Set3:*/
	matrix rowF[1,2]=-2266	/*	Females*/
	
	matrix row2=J(1,18,-6666)
	matrix row2[1,1]=-111		/*	Weight*/
	matrix row2[1,2]=-222		/*	-for-*/
	matrix row2[1,3]=-333		/*	age*/
	
	matrix row11=J(1,18,-6666)
	matrix row11[1,1]=-1111		/*	Height*/
	matrix row11[1,3]=-2222		/*	-for-*/
	matrix row11[1,4]=-3333		/*	age*/
		
	matrix row29=J(1,18,-6666)
	matrix row29[1,1]=-11111	/*	BMI*/
	matrix row29[1,2]=-22222	/*	-for-*/
	matrix row29[1,3]=-33333	/*	age*/
	
	matrix row47=J(1,18,-7777)

	
	matrix Comb=row1 \ row2 \ T \ _zwfa \ row11 \ T \ _zhfa \ row29 \ T \ _zbfa \ row47
	matrix Comb1=rowM \ row2 \ T \ _zwfa1 \ row11 \ T \ _zhfa1 \ row29 \ T \ _zbfa1 \ row47
	matrix Comb2=rowF \ row2 \ T \ _zwfa2 \ row11 \ T \ _zhfa2 \ row29 \T  \ _zbfa2 \ row47
	matrix CombA=Comb \ Comb1 \ Comb2 
	matrix __CombA=first, CombA 
	svmat __CombA
	} 
	/* quietly ends*/
	
	qui gen str3 __xCombA1 = string(__CombA1)
	qui replace __xCombA1 = " " if __xCombA1 == "-6666"
	qui replace __xCombA1 = " " if __xCombA1 == "47" | __xCombA1 == "94" | __xCombA1 == "141"
	
	foreach i in 1 2 11 29 47 48 49 58 76 95 96 97 105 123 {
		qui replace __xCombA1=" " if __xCombA1=="`i'" 
	}
	foreach i in 3 12 30 50 59 77 97 106 124 {
		qui replace __xCombA1="Age" if __xCombA1=="`i'" 
	}
	foreach i in 4 51 98 {
		qui replace __xCombA1="(5-10)" if __xCombA1=="`i'" 
	}
	foreach i in 13 31 60 78 107 125 {
		qui replace __xCombA1="(5-19)" if __xCombA1=="`i'" 
	}
	foreach i in 5 14 32  52 61 79 99 108 126  {
		qui replace __xCombA1="5" if __xCombA1=="`i'" 
	}
	foreach i in 6 15 33  53 62 80 100 109 127  {
		qui replace __xCombA1="6" if __xCombA1=="`i'" 
	}
	foreach i in 7 16 34  54 63 81 101 110 128  {
		qui replace __xCombA1="7" if __xCombA1=="`i'" 
	}
	foreach i in 8 17 35  55 64 82 102 111 129  {
		qui replace __xCombA1="8" if __xCombA1=="`i'" 
	}
	foreach i in 9 18 36  56 65 83 103 112 130  {
		qui replace __xCombA1="9" if __xCombA1=="`i'" 
	}
	foreach i in 10 19 37  57 66 84 104 113 131  {
		qui replace __xCombA1="10" if __xCombA1=="`i'" 
	}
	foreach i in 20 38 58 67 85 107 114 132  {
		qui replace __xCombA1="11" if __xCombA1=="`i'" 
	}
	foreach i in 21 39 59 68 86 108 115 133  {
		qui replace __xCombA1="12" if __xCombA1=="`i'" 
	}
	foreach i in 22 40 60 69 87 109 116 134  {
		qui replace __xCombA1="13" if __xCombA1=="`i'" 
	}
	foreach i in 23 41 61 70 88 110 117 135  {
		qui replace __xCombA1="14" if __xCombA1=="`i'" 
	}
	foreach i in 24 42 62 71 89 111 118 136  {
		qui replace __xCombA1="15" if __xCombA1=="`i'" 
	}
	foreach i in 25 43 63 72 90 112 119 137  {
		qui replace __xCombA1="16" if __xCombA1=="`i'" 
	}
	foreach i in 26 44 64 73 91 113 120 138  {
		qui replace __xCombA1="17" if __xCombA1=="`i'" 
	}
	foreach i in 27 45 65 74 92 114 121 139  {
		qui replace __xCombA1="18" if __xCombA1=="`i'" 
	}
	foreach i in 28 46 66 75 93 115 122 140  {
		qui replace __xCombA1="19" if __xCombA1=="`i'" 
	}
	forval i = 2/19 {
		qui gen str10 __xCombA`i' = string(__CombA`i')
		qui replace __xCombA`i' = "N" if __xCombA`i' == "-1000"
		qui replace __xCombA`i' = "%<-3SD" if __xCombA`i' == "-2000"
		qui replace __xCombA`i' = "%<-2SD" if __xCombA`i' == "-3000"
		qui replace __xCombA`i' = "Mean" if __xCombA`i' == "-4000"
		qui replace __xCombA`i' = "SD" if __xCombA`i' == "-5000"
		qui replace __xCombA`i' = "%>+1SD" if __xCombA`i' == "-6000"
		qui replace __xCombA`i' = "%>+2SD" if __xCombA`i' == "-7000"
		qui replace __xCombA`i' = "%>+3SD" if __xCombA`i' == "-8000"
		qui replace __xCombA`i' = "95%" if __xCombA`i' == "-9000"
		qui replace __xCombA`i' = "C.I." if __xCombA`i' == "-10000"
		qui replace __xCombA`i' = "Set 1:" if __xCombA`i' == "-11"
		qui replace __xCombA`i' = "Sexes" if __xCombA`i' == "-22"
		qui replace __xCombA`i' = "combined" if __xCombA`i' == "-33"
		qui replace __xCombA`i' = "Set 2:" if __xCombA`i' == "-1155"
		qui replace __xCombA`i' = "Males" if __xCombA`i' == "-2255"
		qui replace __xCombA`i' = "Set 3:" if __xCombA`i' == "-1166"
		qui replace __xCombA`i' = "Females" if __xCombA`i' == "-2266"
		qui replace __xCombA`i' = "Weight" if __xCombA`i' == "-111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-222"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-333"
		qui replace __xCombA`i' = "Height" if __xCombA`i' == "-1111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-2222"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-3333"
		qui replace __xCombA`i' = "BMI" if __xCombA`i' == "-11111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-22222"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-33333"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-6666"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-7777"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-69"
	}

	/***	Writing out to the worksheet **/
	local string "xx\yy_prev.xls"
	local i=`datalib'
	local j=`datalab'
	global outp: subinstr local string "xx" "`i'"
	global outprev: subinstr global outp "yy" "`j'"
	outsheet __xCombA1-__xCombA19 using "$outprev" in 1/141, nonames nolabel replace
	
	di " "
	di "Note:	Prevalences are written to"
	di "		$outprev"
	
	/* Cleaning up after Phase II*/
	drop _*
	tempvar clear
end
exit
