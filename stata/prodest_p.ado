*! version 1.0.1 27Sep2016
*! version 1.0.2 05Jun2017 major changes in the code and how the whole routine works, added exponential and parameters options
*! authors: Gabriele Rovigatti, University of Chicago Booth, Chicago, IL & EIEF, Rome, Italy. mailto: gabriele.rovigatti@gmail.com
*!          Vincenzo Mollisi, Bolzano University, Bolzano, Italy & Tor Vergata University, Rome, Italy. mailto: vincenzo.mollisi@gmail.com

/***************************************************************************
** Stata program for prodest Postestimation
**
** Programmed by: Gabriele Rovigatti

**************************************************************************/

cap program drop prodest_p
program define prodest_p, sortpreserve eclass

	version 10.0

	syntax [anything] [if] [in] [, 			   ///
		RESIDuals 				   /// 
		EXPonential 			   ///
		PARameters					///
		]

	marksample touse 	// this is not e(sample)
	
	tempvar esample
	qui gen byte `esample' = e(sample)
	
	loc varlist `anything'
	
	loc mod = "`e(PFtype)'"
	
	/* check for options in launching command */
	if ( "`residuals'" == ""  & "`exponential'" == "" & "`parameters'" == ""){
		di as error "You must specify RESIDuals, EXPonential or PARameters"
		exit 198
	}
	/* check for correct usage of options */
	if ( ("`residuals'" != "" | "`exponential'" != "") & "`parameters'" != ""){
		di as error "the 'parameters' option cannot be used with other options"
		exit 198
	}
	
	if "`mod'" == "Cobb-Douglas"{ /* PART I: COBB-DOUGLAS */
		if ("`residuals'" != "" | "`exponential'" != "") {
			tempname beta
			mat `beta' = e(b)
			tempvar rhs 
			mat score double `rhs' = `beta'
			loc lhs `e(depvar)'
			if "`exponential'" != ""{
				qui gen `varlist' = exp(`lhs' - `rhs') `if'
			}
			else{
				qui gen `varlist' = `lhs' - `rhs' `if'
			}
		}
		else{ /* 'parameters' with cobb-douglas PF yields the results' table */
			di _coef_table, level($S_level)
		}
	}
	else { /* PART II: TRANSLOG */
		loc free = "`e(free)'"
		loc state = "`e(state)'"
		loc controls = "`e(controls)'"
		loc transvars `free' `state' `controls'
		loc translogNum: word count `transvars'
		
		tempname beta
		mat `beta' = e(b) // extract the estimated betas
		
		loc n = 1 // regenerate the variables used in the routine in order to fit the values
		foreach x of local transvars{
			tempvar var_`n' betavar_`n' 
			qui g `betavar_`n'' = `beta'[1,`n'] * `x'
			qui g `var_`n'' = `x'
			loc fit `fit' -`betavar_`n''
			loc ++n
		}
		forv i = 1/`translogNum'{
			forv j = `i'/`translogNum'{ /* `i' */
				tempvar var_`i'`j' betavar_`i'`j'
				cap g `betavar_`i'`j'' = `beta'[1,`n'] * (`var_`i'' * `var_`j'')
				cap g `var_`i'`j'' = (`var_`i'' * `var_`j'')
				loc ++n
			}
		}
		if "`exponential'" != "" {
			qui g `varlist' = exp(`e(depvar)' `fit') `if' // here generate the predicted residuals -- exponential
		}
		else if "`residuals'" != ""{
			qui g `varlist' = `e(depvar)' `fit' `if' // here generate the predicted residuals
		}
		else{ /* in case of 'parameters' option */
			loc freenum: word count `free'
			loc statenum: word count `state'
			loc totnum:  word count `free' `state' 
			forv i = 1/`totnum'{
				forv j = 1/`totnum'{
					if `i' != `j'{ /* generate the cross variables part only  */
						cap confirm variable `betavar_`i'`j''
						if !_rc{
							loc remainder `remainder' + (`betavar_`i'`j''/`var_`i'')
						}
					}
				}
				tempvar betafit_`i' // the parameter for translog is defined as beta_Wtranslog = beta_w + 2*beta_ww * W + beta_wx * X 
				qui gen `betafit_`i'' = `beta'[1,`i'] + 2*(`betavar_`i'`i''/`var_`i'') `remainder' // here we use the previously generated variables and weight them by the ith variable
				qui su `betafit_`i'', meanonly
				loc beta_`i': di %6.3f `r(mean)'
				loc remainder ""
			}
		di _n _n
		di as text "{hline 75}"
		di as text "Translog elasticity estimates" _continue
		di _col(49) "prodest postestimation"
		di as text "{hline 75}"
		di as text "Elasticity Parameter" _continue
		di _col(49) "Value"
		di as text "{hline 75}"
		loc i = 1
		foreach var of varlist `free' `state'{
			di as text "beta_`var'" _continue
			di _col(49) "`beta_`i''"
			loc ++i
		}
		di as text "{hline 75}"
		}
	}

end
