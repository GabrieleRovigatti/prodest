*! version 1.0.1 10May2020
*! author: Gabriele Rovigatti, Bank of Italy, Rome, Italy. mailto: gabriele.rovigatti@gmail.com  |  gabriele.rovigatti@bancaditalia.it
/***************************************************************************
** Stata program for Markup Estimation - prodest postestimation
** Programmed by:   Gabriele Rovigatti
**************************************************************************/

cap program drop prodest_p_m
program define prodest_p_m, sortpreserve eclass

	version 10.0

	syntax [anything] [if] [in] [, 			   ///
		MARKups						///
		INPUTvar(varlist numeric min=1 max=1)		///
		CORRected					///
		REPetition(integer 1)			///
		]

	marksample touse 	// this is not e(sample)
	
	tempvar esample
	qui gen byte `esample' = e(sample)
	
	loc varlist `anything'
	
	loc mod = "`e(PFtype)'"
	loc fsres = "`e(FSres)'"
	loc free = "`e(free)'"
	tempvar val
	if "`e(model)'" == "grossoutput"{
		loc proxy = "`e(proxy)'"
		loc free "`free' `proxy'"
		qui g `val' = log( exp(`e(depvar)') - exp(`proxy') )  // in case of gross output, generate the measure of value added as the difference between gross output and material input costs
	}
	else{
		qui g `val' = `e(depvar)'
	}
	loc state = "`e(state)'"
	
	************************************************************************************************
			********************* PART 1: MARKUPS ESTIMATION ************************
	************************************************************************************************s
	/* check for correct usage of options */
	if !mi("`corrected'") & mi("`fsres'"){ // correction à la DLW only available with first-stage residuals
		di as error "Markup correction requires launching prodest with 'fsresiduals(<fsvar>)' option"
		exit 198	
	}
	if !`:list inputvar in free'{ // check whether the input variable specified is in the list of free variables
		di as error "<inputvar> should be either a free or a proxy variables used in the estimation"
		exit 198
	}
	cap confirm var `varlist', exact // if the variable already exists
	if !_rc & `repetition' == 1{
		di as error "`varlist' already exists"
		exit 198
	}
	if mi("`varlist'"){ // check the outcome variable: if it is missing, use a pre-specified _mkup[number] variable
		loc c = 0
		while (1){
			loc ++c
			loc varlist _mkup`c'
			cap confirm var `varlist'
			if (_rc != 0) continue, break
		} 
		di as error "You should specify <newvarname> to store the estimated markups. They will be stored in `varlist' now"
	}
	
	********* ROUTINE START **************
	tempvar theta alpha
	/* generate the input share, either "raw" or corrected by the first-stage residuals as suggested by DLW */
	*loc lhs `e(depvar)' // this is the output - in logs
	qui g `alpha' = exp(`inputvar') / exp(`val') // share is input cost / value added
	if !mi("`corrected'"){
		qui replace `alpha' = `alpha' * exp(`fsres')
	}
	/* Generate the elasticity parameter - either the estimated beta (Cobb-Douglas) or a function of it (Translog) */ 
	if "`mod'" == "Cobb-Douglas"{ /* PART I: COBB-DOUGLAS */
		qui g `theta' = _b[`inputvar'] // 
	}
	else { /* PART II: TRANSLOG */
		tempname beta
		mat `beta' = e(b) // extract the estimated betas

		loc controls = "`e(controls)'"
		loc transvars `free' `state' `controls'
		loc translogNum: word count `transvars'
		
		loc n = 1 // regenerate the variables used in the routine in order to fit the values
		foreach x of local transvars{
			tempvar var_`n' betavar_`n' 
			qui g `var_`n'' = `x'
			loc ++n
		}
		forv i = 1/`translogNum'{
			forv j = `i'/`translogNum'{ 
				tempvar var_`i'`j' beta_`i'`j'
				cap g `var_`i'`j'' = (`var_`i'' * `var_`j'')
				cap g `beta_`i'`j'' = `beta'[1,`n']
				loc ++n
			}
		}
		loc varnum:  word count `free' `state' 
		loc inputpos: list posof "`inputvar'" in free // this is the number of the input variable within the free vars
		forv j = 1/`varnum'{
			if `inputpos' != `j'{ /* generate the cross variables part only  */
				cap confirm variable `beta_`inputpos'`j''
				if !_rc{
					loc remainder `remainder' + (`beta_`inputpos'`j'' * `var_`j'') 
				}
				else{
					loc remainder `remainder' + (`beta_`j'`inputpos'' * `var_`j'')
				}
			}
		}
		qui gen `theta' = `beta'[1,`inputpos'] + (2 * `beta_`inputpos'`inputpos'' * `var_`inputpos'') `remainder' // the elasticity for translog is defined as beta_Wtranslog = beta_w + 2*beta_ww * W + beta_wx * X, and here we use the previously generated variables and weight them by the ith variable
	}
	/* Compute the Markups */
	if `repetition' > 1{
		tempvar _foo_
		g `_foo_' = `theta' / `alpha' `if'
		replace `varlist' = `_foo_' `if'
	}
	else{
		g `varlist' = `theta' / `alpha' `if' // compute the markups, and save the relative variable
	}
end
