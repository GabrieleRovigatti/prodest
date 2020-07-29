*! version 1.0.1 10May2020
*! author: Gabriele Rovigatti, Bank of Italy, Rome, Italy. mailto: gabriele.rovigatti@gmail.com  |  gabriele.rovigatti@bancaditalia.it
/***************************************************************************
** Stata program for Markup Estimation 
** Programmed by:   Gabriele Rovigatti
**************************************************************************/
/* 
Command structure: 
- MICRO
	- 1) DLW
		Estimation Type:
			|- PRODEST (CUSTOM VERSION)
			|- SAVED PRODEST RESULTS
		- PREDICT MARKUPS
- MACRO
	- 2.1) HALL (1988,2018) 
		Class:
			|- GROSS OUTPUT
		Options and Postestimation:
			|- TIME-VARYING
			|- MEASUREMENT ERROR CORRECTION
	- 2.2) ROEGER (1995)
	
*/
// ********************* //

capture program drop markupest
program define markupest, sortpreserve eclass byable(recall)
	version 10.0
	syntax name [if] [in] [, METhod(name) replace SAVEBeta VERBose id(varlist min=1 max=1) t(varlist min=1 max=1) output(varlist numeric min=1 max=1) free(varlist numeric min=1) /*
		*/ state(varlist numeric min=1 max=1) proxy(varlist numeric min=1) VALueadded PRODESTOPTions(string) CORRected ESTimates(name) INPUTvar(varlist numeric min=1) /*
		*/ GO(varlist numeric min=1 max=1) deltago(varlist min=1 max=1) pgo(varlist min=1 max=1) INputs(varlist numeric min=1) INSTRuments(varlist numeric min=1) /*
		*/ DELTAvars(varlist numeric min=1) PRICEvars(varlist numeric min=1) TIMEVarying HGraph pmethod(name) W(varlist numeric min=1 max=1)] 
	
	loc vv: di "version " string(min(max(10,c(stata_version)),14.0)) ", missing:" // we want a version 10 min, while a version 14.0 max (no version 14.2)
	loc markupvar `namelist'
	marksample touse
	markout `touse' `output' `free' `state' `proxy' `go' `va' `compensation' `ii' `rt' `k' `w' `instruments' 
	
	**************** COMMAND CHECKS *********************
	if !inlist("`method'", "dlw", "hall", "roeger"){ // if the users specifies a wrong method
		di as error "Allowed methods are dlw, hall and roeger"
		exit 198
	}
	/* Indicate either Micro or Macro methods */
	loc micro = 0
	if "`method'" == "dlw" | !mi("`output'") | !mi("`free'") | !mi("`state'") | !mi("`proxy'") | !mi("`estimates'"){
		loc micro = 1
	}
	/* method definition in case of missing */
	if mi("`method'") & `micro' == 1{ // in case of micro the only available method is DLW
		di as error "'method()' is missing. The default 'dlw' (micro) is used."
		loc method dlw 
	}
	else if mi("`method'") & `micro' == 0{ // in case of macro, Roeger is the default
		di as error "'method()' is missing. The default 'roeger' (macro) is used."
		loc method roeger
	}
	/* check whether data is xtset. In case it is not, xtset it with the specified variables */
	if (mi("`id'") | mi("`t'")) {
	  capture xtset
	  if (_rc != 0) { // if not xtset, the user must specify both panel and time variables
		 di as error "You must either xtset your data or specify both id() and t()"
		 error _rc
	  }
	  else {
		 loc id = r(panelvar)
		 loc t = r(timevar)
	  }
	}
	else {
	  qui xtset `id' `t'
	}
	
	loc byRep = _byindex() // generate an indicator for the by repetition. If it is not the first, we should take care of that
	
	/* check whether the variable already exists */
	cap confirm variable `markupvar', exact // check whether the markup variable exists already
	if _rc == 0 & `byRep' == 1{
		if !mi("`replace'"){ // if the user specified "replace", drop the current variable
			drop `markupvar'
		}
		else{
			di as error "`markupvar' already exists"
			exit 198
		}
	}
	else if _rc != 0 & `byRep' == 1 & !mi("`replace'"){ // in case the user specifies "replace", but there is no variable
		di as error "`markupvar' does not exist. No replace"
	}
	/* ROUTINE START  */
	if `micro' == 1{
		if mi("`inputvar'"){
			di as error "You should specify the <inputvar> whose elasticity would be used to compute markups"
			exit 198
		}
		************************************************************************************
		**************** MICRO PART: PRODEST + OUTPUT SHARE CALCULATION ********************
		************************************************************************************
		/* 1) DLW */
		if !mi("`output'") | !mi("`free'") | !mi("`state'") | !mi("`proxy'"){ // if the user specifies at least one of the variables, this is the "micro" run
			if mi("`output'"){ // check whether there are missing "parts" in the command launched
				di as error "Specify an <output> variable "
				exit 198
			}
			else if mi("`free'"){
				di as error "Specify at least one <free> variable "
				exit 198
			}
			else if mi("`state'"){
				di as error "Specify at least one <state> variable "
				exit 198
			}
			else if mi("`proxy'"){
				di as error "Specify at least one <proxy> variable "
				exit 198
			}
			if mi("`pmethod'"){ // the default method is levisohn-petrin
				loc pmethod lp
			}
			if !mi("`corrected'") & !regexm("fsres", "`prodestoptions'") & !inlist("`method'", "wrdg", "mr", "rob") { // if the user specifies "corrected markups", but not the first-stage residuals, run them with a tmp var
				tempvar fsresid
				loc fscorr fsresiduals(`fsresid')
			}
			/* run prodest with the specified variables + options */
			prodest_m `output' if `touse' == 1, free(`free') state(`state') proxy(`proxy') met("`pmethod'") `valueadded' `fscorr' `prodestoptions' 
			if !mi("`verbose'"){
				_coef_table
			}
			if !mi("`savebeta'"){
				mat beta = e(b)
				loc betanum = colsof(beta)
				loc betanames: colnames beta
				forv v = 1/`betanum'{
					loc betaname `: word `v' of `betanames''
					qui g _b`betaname' = beta[1, `v']
				}
			}
			/* run the PREDICT */
			if !mi("`e(FSres)'"){ // if prodest is run with first-stage residuals, correct for them - otherwise go with the uncorrected 
				qui predict `markupvar' if `touse' == 1, markups inputvar(`inputvar') corrected rep(`byRep') 
			}
			else{
				qui predict `markupvar' if `touse' == 1, markups inputvar(`inputvar') rep(`byRep') 
			}
		}
		else if !mi("`estimates'") { // estimates stored
			estimates restore `estimates'	
			loc free `e(free)'
			loc state `e(state)'
			loc proxy `e(proxy)'
			loc output `e(depvar)'			
			/* run the PREDICT */
			if !mi("`e(FSres)'"){ // if prodest is run with first-stage residuals, correct for them - otherwise go with the uncorrected 
				qui predict `markupvar' if `touse' == 1, markups inputvar(`inputvar') corrected rep(`byRep')
			}
			else{
				qui predict `markupvar' if `touse' == 1, markups inputvar(`inputvar') rep(`byRep') 
			}
		}
	}
	else{
		******************************************************************************************
		**************** MACRO PART: HALL (1988), ROEGER (1995) and HALL (2018) ******************
		******************************************************************************************
		/* Listing macro-specific errors */
		if mi("`go'"){
			di as error "You must specify <GO> variable"
			exit 198
		}
		loc nInputs: word count `inputs' // check that there are at least capital and labor inputs specified
		if mi("`inputs'") | `nInputs' < 2{
			di as error "You must specify at least two <INputs> (K and L)"
			exit 198
		}
		if "`method'" == "hall" {
			**********************************************************************
			/* 2.1) Hall method */
			**********************************************************************
			*** Sanity checks ***
			if mi("`instruments'"){
				di as error "hall method requires at least one <INSTRument>"
				exit 198
			}
			if !mi("`valueadded'") & `nInputs' < 3{ // with Value Added models, there must be at least 3 inputs (Labor, Capital, and Materials)
				di as error "With VA models, you must specify at least three <INputs> (K, L, and M)"
				exit 198
			}
			if mi("`pricevars'") & mi("`deltavars'"){ // input 
				di as error "You should specify either <DELTAvars> or <PRICEvars>"
				exit 198
			}
			loc nDeltaPrice: word count `deltavars' `pricevars' // check that either deltavars of pricevars are consistent with the input variables by checking that the sum of the number is consistent
			if `nDeltaPrice' != `nInputs'{
				di as error "The number of <DELTAvar> of <PRICEvar> specified must equal the <INputs>"
			}
			// Left-hand side: sum of shares + changes in input
			if mi("`deltavars'"){ // generate deltavars --> requires the "pricevars" - i.e., the changes in prices
				loc i = 0
				foreach var of varlist `inputs'{
					tempvar d`var' delta`var'
					loc i = `i' + 1
					loc pvar: word `i' of `pricevars'
					qui g `d`var'' = (`var' - l.`var') / l.`var' `if' // generate the changes in the "gross" variable
					qui g `delta`var'' = `dvar' - `pvar' // subtract the changes in prices to obtain the changes in input
					loc deltavars `deltavars' `delta`var''
				}
			}
			// generate the sum of weighted inputs
			loc i = 0 //counter for inputs
			tempvar lhs
			qui g `lhs' = 0
			foreach dvar of varlist `deltavars'{
				tempvar _alpha`var' 
				loc i = `i' + 1
				loc invar: word `i' of `inputs' // this is the level variable relative to the delta we are considering
				qui g `_alpha`var'' = `invar' / `go' // generate the alphas of input variables - that is, the ratio between the level of inputvar cost, and gross output
				qui replace `lhs' = `lhs' + `_alpha`var'' * `dvar' // sum of weighted input variables changes
			}
			// Right-hand side: changes in output
			tempvar rhs
			if !mi("`deltago'"){
				qui g `rhs' = `deltago'
			}
			else{ // generate the changes in sectoral output / value added
				tempvar dGO deltaGO
				qui g `dGO' = (`go' - l.`go') / l.`go' // this is the change in Gross Output = P * Y
				qui g `deltago' = `dGO' - `pgo' // here we generate delta Y 
				qui g `rhs' = `deltago'
			}
			// Time-varying markups
			if !mi("`timevarying'"){
				tempvar twdeltaGO tweights
				qui su `t' `if'
				loc avgY = `r(min)' + int( (`r(max)' - `r(min)') / 2 )
				g `tweights' = `t' - `avgY' // linear time trend
				qui g `twdeltaGO' = `tweights' * `deltago'
				foreach var of varlist `instruments'{
					tempvar tw`var'
					qui g `tw`var'' = `tweights' * `var'
					loc twinstruments `twinstruments' `tw`var''
				}
			}
			// Actual estimation
			qui ivreg2 `lhs' (`rhs' `twdeltaGO' = `instruments' `twinstruments')  if `touse' == 1, nocons
			if `byRep' == 1{ // when using "by()", discriminate between the first and the following rounds
				if mi("`timevarying'"){
					qui g `markupvar' = 1 / _b[`rhs'] if `touse' == 1
				}
				else{
					qui g _psi_ = -_b[`twdeltaGO'] if `touse' == 1
					qui g `markupvar' = 1 / _b[`rhs'] - (_b[`twdeltaGO'] * `tweights') if `touse' == 1
				}
			}
			else{ // 
				tempvar _tmp_ 
				if mi("`timevarying'"){
					qui g `_tmp_' = 1 / _b[`rhs'] if `touse' == 1
				}
				else{
					qui replace _psi_ = -_b[`twdeltaGO'] if `touse' == 1
					qui g `_tmp_' = 1 / _b[`rhs'] - (_b[`twdeltaGO'] * `tweights') if `touse' == 1
				}
				qui replace `markupvar' = `_tmp_' if `touse' == 1
			}
			// post-estimation cleaning
			if !mi("`hgraph'") & _bylastcall() == 1{ // plot the resulting, corrected distribution - if it is the last "by" call - or the unique one
				preserve
					tempvar densT densWerror
					collapse (mean) `markupvar', by(`id')
					qui drop if `markupvar' < 0 // way to ensure that initial sigma > 0 --> which ensures that the estimation is feasible
					qui g muMinus1 = `markupvar' - 1
					cap qui mata: opt_hall()
					clear
					set obs 1000000
					qui g `densT' = 1 + exp(rnormal(delta, sigma)) // 
					qui g `densWerror' = rnormal(meanM, sdTotM ) // 
					qui replace `densT' = 0 if _n == 1
					tw (kdensity `densT' if `densT' >= 0 & `densT' < 3, lw(thick) lc(maroon)) (kdensity `densWerror' if `densWerror' >= 0 /*
						*/ & `densWerror' < 3, lw(thick) lc(ebblue) ), xtitle("Ratio of price to marginal cost {&mu}") xlab(0(0.25)3) /*
						*/ legend(order(1 "{&mu} w/o sampling error" 2 "Estimated {&mu}")) ytitle("Probability Density")
				restore
			}
		}
		else if "`method'" == "roeger"{
			**********************************************************************
			/* 2.2) Roeger method */
			**********************************************************************
			if mi("`go'"){ // check whether there are missing "parts" in the command launched
				di as error "Specify the <GO> variable "
				exit 198
			}
			if !mi("`w'"){
				loc weights "[aw=`W']"			
			}
			loc nInputs: word count `inputs' // this is the total number of inputs
			
			tempvar dGO
			qui g `dGO' = (`go' - l.`go') / l.`go' // changes in gross output
			
			/* generate the deltavars for Roeger: SUBSTITUTE THE PREVIOUS CODE */
			foreach var of varlist `inputs'{
				tempvar d`var'
				qui g `d`var'' = (`var' - l.`var') / l.`var'
				loc ddvars `ddvars' `d`var''
			}
				
			// generate the sum of weighted inputs
			loc i = 0 //counter for inputs
			tempvar lhs allalphas 
			qui g `lhs' = `dGO'
			qui g `allalphas' = 0
			*foreach dvar of varlist `deltavars'{
			foreach dvar of varlist `ddvars'{
				tempvar _alpha`var' 
				loc i = `i' + 1
				loc invar: word `i' of `inputs' // this is the level variable relative to the delta we are considering
				if `i' < `nInputs'{ // generate alphas for each input but capital, which has the "residual" share
					qui g `_alpha`var'' = `invar' / `go' // generate the alphas of input variables - that is, the ratio between the level of inputvar cost, and gross output
					qui replace `allalphas' = `allalphas' + `_alpha`var'' 
				}
				else{ // it should NOT change anything if the user has complete data
					qui g `_alpha`var'' = 1 - `allalphas' if `touse' == 1 // this is the residual alpha for capital
				}
				qui replace `lhs' = `lhs' - `_alpha`var'' * `dvar' // sum of weighted input variables changes
			}
			// Right-hand side: changes in output
			tempvar rhs deltaRtK
			loc capitalvar: word `nInputs' of `inputs'
			qui g `deltaRtK' = (`capitalvar' - l.`capitalvar') / l.`capitalvar' `if' // this is the gross change in capital
			qui g `rhs' = `dGO' - `deltaRtK'
			
			// Actual estimation
			qui reg `lhs' `rhs' `weights' if `touse' == 1, nocons
			if `byRep' == 1{ // when using "by()", discriminate between the first and the following rounds
				qui g `markupvar' = 1 / (1 - _b[`rhs']) if `touse' == 1
			}
			else{ // 
				tempvar _tmp_ 
				qui g `_tmp_' = 1 / (1 - _b[`rhs']) if `touse' == 1
				qui replace `markupvar' = `_tmp_' if `touse' == 1
			}
		}
	}
	
	/* return the locals in e() */
	eret clear 
	
	eret loc cmd "markupest"
	eret loc markupvar "`markupvar'"
	eret loc method "`method'"
	eret loc id "`id'"
	eret loc t "`t'"
	if `micro' == 1 {
		eret loc markuptype "micro"
		eret loc inputvar "`inputvar'"
		if !mi("`estimates'"){
			eret loc output "`output'"
			eret loc free "`free'"
			eret loc state "`state'"
			eret loc proxy "`proxy'"
			eret loc PFest_method "`pmethod'"
		}
		else{
			eret loc estimate_name "`estimates'"
		}
	}
	else{
		eret loc markuptype "macro"
		eret loc inputs "`inputs'"
		if !mi("`instruments'"){
			eret loc instruments "`instruments'"
		}
		if !mi("`deltavars'"){
			eret loc deltavars "`deltavars'"
			eret loc deltago "`deltago'"
		}
		else if !mi("`pricevars'"){
			eret loc pricevars "`pricevars'"
			eret loc pgo "`pgo'"
		}
		if !mi("`timevarying'"){
			eret loc hall_mkuptype "time-varying"
		}
	}
	if !mi("`prodestoptions'"){
		eret loc prodestOPTS "`prodestoptions'"
	}
	if !mi("`corrected'"){
		eret loc corrected "corrected"
	}
	
end program


/*---------------------------------------------------------------------*/
capture mata mata drop hallsolver()
capture mata mata drop opt_hall()


mata:
/*------------ MATA ROUTINE FOR HALL (2018) POSTESTIMATION ------------*/

	void hallsolver(todo, p, M, lnf, S, H)
	{
		 gamma = p[1]
		 delta = p[2]
		 sigma = p[3]
		 
		 M1 = M[1]
		 M2 = M[2]
		 M3 = M[3]
		 
		 lnf = ( exp(delta + 0.5 * (sigma^2) ) - M1)^2 \   
			   ( (gamma^2) + exp( (2 * delta) + (2 * sigma^2) ) - M2 )^2 \  
			   ( exp( (3 * delta) + (9/2) * sigma^2 )  + 3 * gamma^2 * exp( delta + (1/2) * sigma^2 ) - M3 )^2 
	}

/*---------------------------------------------------------------------*/

	void opt_hall()
	{
		st_view(MUminus1=., ., "muMinus1")
		MUminus1_sq = MUminus1:^2
		MUminus1_cb= MUminus1:^3
		M = mean(MUminus1), mean(MUminus1_sq), mean(MUminus1_cb)
		S = optimize_init()
		optimize_init_argument(S, 1, M)
		optimize_init_evaluator(S, &hallsolver())
		optimize_init_evaluatortype(S, "v0")
		optimize_init_params(S, (1, 1, 1) )
		optimize_init_which(S, "min" )
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_ptol(S, 1e-16)
		optimize_init_conv_vtol(S, 1e-16)
		p = optimize(S)
		// these are the three main elements
		gamma = p[1]
		delta = p[2]
		sigma = p[3]
		// compute means and variances - etas, nus and markups
		meanV = exp( delta + (sigma^2 / 2) ) // E[ exp( log(v) ) ] = exp( mu + sigma^2 / 2 )
		varV = ( exp(sigma^2 ) - 1 ) * exp(2 * delta + sigma^2) // Var[ exp( log(v) ) ] = ( exp(sigma^2 ) - 1 ) exp(2 * mu + sigma^2)
		sdV = (varV)^(1/2)
		meanM = 1 + meanV
		sdM = sdV // the standard deviation of the markup is the standard deviation of nu, i.e., the estimated markups net of the disturbance
		sdTotM = (varV + gamma^2)^(1/2) // this is the standard deviation of the estimated markups
		
		st_numscalar("delta", delta)
		st_numscalar("sigma", sigma)
		st_numscalar("meanM", meanM)
		st_numscalar("sdTotM", sdTotM)
	}

/*---------------------------------------------------------------------*/
end



