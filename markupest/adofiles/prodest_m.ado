**! version 1.0.1 10May2020
*! author Gabriele Rovigatti, Bank of Italy, Rome, Italy. mailto: gabriele.rovigatti@gmail.com  |  gabriele.rovigatti@bancaditalia.it

/***************************************************************************
** Stata program for Markup Estimation - custom PRODEST version
** Programmed by:   Gabriele Rovigatti
**************************************************************************/

capture program drop prodest_m
program define prodest_m, sortpreserve eclass
	version 10.0
	syntax varlist(numeric min=1 max=1) [if] [in], free(varlist numeric min=1) /*
		  */ proxy(varlist numeric min=1 max=2) state(varlist numeric min=1) /*
		  */ [control(varlist min=1) ENDOgenous(varlist min=1) id(varlist min=1 max=1) t(varlist min=1 max=1) reps(integer 5) /*
		  */ VAlueadded Level(int ${S_level}) OPTimizer(namelist min=1 max=3) MAXiter(integer 10000) /* OPTIMIZER: THERE IS THE POSSIBILITY TO WRITE technique(nr 100 nm 1000) with different optimizers after the number of iterations. 
		  */ poly(integer 3) METhod(name min=1 max=1) lags(integer 999) TOLerance(real 0.00001) gmm /*
		  */ ATTrition ACF INIT(string) FSRESiduals(name min=1 max=1) EVALuator(string) TRANSlog OVERidentification]  
	
	loc vv: di "version " string(min(max(10,c(stata_version)),14.0)) ", missing:" // we want a version 10 min, while a version 14.0 max (no version 14.2)
	loc depvar `varlist'
	marksample touse
	markout `touse' `free' `proxy' `state' `control' `endogenous' `id' `t' 
	
	/// checks for same variables in state - free - proxy - controls
	loc chk1: list free & state
	loc chk2: list free & control
	loc chk3: list free & proxy
	loc chk4: list state & control
	loc chk5: list state & proxy
	loc chk6: list control & proxy
	forval i = 1/6{
		if "`chk`i''" != ""{
			di as error "Same variables in free, state, control or proxy."
			exit 198
		}
	}
	/// check whether there are more than one state AND more than one proxy
	loc pnum: word count `proxy'
	loc snum: word count `state'
	if `pnum' > 1 & `snum' > 1{
		di as error "Cannot specify multiple state AND multiple proxy variables"
		exit 198
	}
	/// check for unavailable choice of models
	if (!inlist("`method'", "op", "lp", "wrdg", "mr", "rob") & !mi("`method'")){
		di as error "Allowed prodest methods are op, lp, wrdg, mr or rob. Default is lp"
		exit 198
	}
	else if ("`method'" == "mr" & c(stata_version) < 14.2){
		di as error "MrEst only available with Stata version 14.2 or higher"
		exit 198
	}
	else if mi("`method'"){
		loc method "lp"
	}
	/// Check for unpractical value of polynomial approximation
	if (`poly' >= 7 | `poly' < 2){
		di as error "Polynomial degree must lie between 2 and 6"
		exit 198
	}
	/// Syntax check: is data xtset? 
	if (mi("`id'") | mi("`t'")) {
	  capture xtset
	  if (_rc != 0) {
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
	/// Check for a valid confidence level 
	if (`level' < 10 | `level' > 99) {
		di as error "confidence level must be between 10 and 99"
		error 198
	}
	/// Value added or gross output? 
	if mi("`valueadded'"){
		loc model "grossoutput"
		loc proxyGO `proxy'
		loc lproxyGO l_gen`proxy'
	}
	else {
		loc model "valueadded"
	}
	/// Number of repetitions reasonable? 
	if (`reps' < 2) & ("`method'" != "wrdg" & "`method'" != "mr"){
		di as error "reps() must be at least 2"
		exit 198
	}
	/// optimizer choice
	if (!mi("`optimizer'") & !inlist("`optimizer'", "nm", "nr", "dfp", "bfgs", "bhhh", "gn")){
		di as error "Allowed optimizers are nm, nr, dfp, bfgs and bhhh or gn for wrdg and mr. Default is nm or gn for wrdg and mr"
		exit 198
	}
	else if inlist("`method'", "wrdg", "mr") & inlist("`optimizer'", "nm", "bhhh"){
		di as error "`optimizer' is not allowed with `method' method. Optimizer switched to gn (default)"
		loc optimizer "gn"
	}
	else if mi("`optimizer'") & !inlist("`method'", "wrdg", "mr"){
		loc optimizer "nm"
	}
	else if mi("`optimizer'"){
		loc optimizer "gn"
	}
	/// number of lags in MR metodology
	if ("`method'" == "mr" & mi("`lags'")){
		di as error "Method 'mr' requires a specification for lags. Switched to 'all possible lags' (default)"
		loc lags = .
	}
	else if ("`method'" == "mr" & (`lags' < 1)) {
		di as error "Minimum lag is 1"
		exit 198
	}
	else if ("`method'" == "mr" & (`lags' == 999)) {
		loc lags = .
	}
	/// check ACF correction in woolridge and mr cases
	if ("`method'" == "wrdg" | "`method'" == "mr") & !mi("`acf'"){
		di as error "`method' does not support ACF correction. Estimation run on baseline model"
		loc acf ""
	}
	/// feasible values of tolerance
	if mi(`tolerance') | `tolerance' > 0.01{
		di as error "maximum value of tolerance is 0.01. Changed to default of value of e-5"
		loc tolerance = 0.00001
	}
	if ("`method'" == "op") {
		loc proxyGO ""
		loc lproxyGO ""
	}
	/// define conv_nrtol for wooldridge and mr - NOT in case of GN
	if ("`method'" == "wrdg" | "`method'" == "mr") & "`optimizer'" != "gn"{
		loc conv_nrtol "conv_nrtol(`tolerance')"
	}
	if ("`method'" == "wrdg" | "`method'" == "mr") & !mi("`init'"){
		loc init = subinstr("`init'",","," ",.)
		foreach var in `free' `state'{
			gettoken val init: init
			loc init_gmm `init_gmm' xb_`var' `val'
		}
		loc init_gmm from(`init_gmm')
	}
	if !mi("`fsresiduals'") & inlist("`method'", "wrdg", "mr", "rob"){
		di as error "fsresiduals is available with OP and LP methods only."
		exit 198
	}
	cap confirm var `fsresiduals'
	if !_rc{
		di as error "`fsresiduals' already exists"
		exit 198
	}
	/// define the evaluator type
	if mi("`evaluator'") & "`optimizer'" == "bhhh"{
		loc evaluator = "gf0"
	}
	else if mi("`evaluator'"){
		loc evaluator = "d0"
	}
	/// check the translog - only meaningful for ACF and Wooldridge methods
	if !mi("`translog'") & mi("`acf'") & "`method'" != "wrdg"{
		di as error "translog is available with Wooldridge or ACF-corrected models only"
		exit 198
	}
	else if !mi("`translog'") & !mi("`acf'"){
		loc transVars `free' `state' `proxyGO'
	}
	/// change the production function type
	if !mi("`translog'"){
		loc PFtype "translog"
	} 
	else{
		loc PFtype "Cobb-Douglas"
	}
	if !mi("`overidentification'") & mi("`acf'") & "`method'" != "wrdg"{
		di as error "overidentification is meaningful for ACF or WRDG methods only. Ignoring the option."
	}
	if !mi("`gmm'") & "`method'" != "wrdg"{
		di as error "gmm works for wrdg only. Ignoring the option."
		loc gmm ""
	}
	if !mi("`acf'") & "`model'" == "grossoutput"{
		di as error "Using ACF correction with GO output does not ensure a correct parameter identification. See ACF (2015)."
	}
	
	if "`method'" == "op" loc strMethod "Olley-Pakes"
	if "`method'" == "lp" loc strMethod "Levinsohn-Petrin"
	if "`method'" == "wrdg" loc strMethod "Wooldridge"
	if "`method'" == "rob" loc strMethod "Wooldridge/Robinson"
	if "`method'" == "mr" loc strMethod "Mollisi-Rovigatti"
	
	loc colnum: word count `free' `state' `control' `proxyGO' 
	
	/// initialize results matrices
	tempname firstb __b __V robV
	mat `__b' = J(`reps',`colnum',.)
	
	/// preserve the data before generating tons of variables
	preserve
		qui su `t' if `touse' == 1
		loc maxDate = `r(max)'
		
		/// directly keep only observations in IF and IN --> SAVE A TEMPORARY FILE 
		tempfile temp
		qui keep if `touse' == 1
		keep `depvar' `free' `state' `proxy' `control' `id' `t' `touse'  `endogenous' 
		/// generate an "exit" dummy variable equal to one for all firms not present in the last period of panel
		tempvar exit
		qui bys `id' (`t'): g `exit' = (_n == _N & `t' < `maxDate')
		qui save `temp', replace
		
		/// define all locals to run the command
		loc toLagVars `free' `state' `control' `proxyGO' 
		
		foreach var of local toLagVars{
			qui g l_gen`var' = l.`var'
			loc laggedVars `laggedVars' l_gen`var'  
		}
		foreach local in free state proxy control{
			loc `local'num: word count ``local''
			foreach var of local `local'{
				loc lag`local'Vars `lag`local'Vars' l_gen`var'
			}
		}
		
		loc instrumentVars `state' `lagfreeVars' `control' `lproxyGO'
		
		/// OP LP polyvars
		loc polyvars `state' `proxy'
		/// ACF requires free variables to be among the polynomial
		if !mi("`acf'"){
			loc polyvars `free' `polyvars' 
		}
		loc varnum: word count `polyvars'
		loc controlnum: word count `control'
		loc tolagnum: word count `toLagVars'
		
		// poly-th degree polynomial
		loc n = 1
		foreach x of local polyvars{
			qui g var_`n' = `x'
			loc interactionvars `interactionvars' var_`n'
			loc ++n
		}

		forv i=1/`varnum'{
			forv j=`i'/`varnum'{
				qui g var_`i'_`j' = var_`i'*var_`j'
				loc interactionvars `interactionvars' var_`i'_`j'
				if `poly' > 2{
					forv z=`j'/`varnum'{
						qui g var_`i'_`j'_`z' = var_`i'*var_`j'*var_`z'
						loc interactionvars `interactionvars' var_`i'_`j'_`z'
						if `poly' > 3{
							forv g = `z'/`varnum'{
								qui g var_`i'_`j'_`z'_`g' = var_`i'*var_`j'*var_`z'*var_`g'
								loc interactionvars `interactionvars' var_`i'_`j'_`z'_`g'
								if `poly' > 4{
									forv v = `g'/`varnum'{
										qui g var_`i'_`j'_`z'_`g'_`v' = var_`i'*var_`j'*var_`z'*var_`g'*var_`v'
										loc interactionvars `interactionvars' var_`i'_`j'_`z'_`g'_`v'
										if `poly' > 5{
											forv s = `v'/`varnum'{
												qui g var_`i'_`j'_`z'_`g'_`v'_`s' = var_`i'*var_`j'*var_`z'*var_`g'*var_`v'*var_`s'
												loc interactionvars `interactionvars' var_`i'_`j'_`z'_`g'_`v'_`s'
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		if inlist("`method'","wrdg", "mr", "rob"){
			/// generate GMM fit - for each state and free variables need to initiate the GMM
			foreach element in `free' `state' `proxyGO'{
				local gmmfit `gmmfit'-{`element'}*`element'
			}
			/// generate lag interactionvars
			foreach var of local interactionvars{
				cap g l_gen`var' = l.`var'
				loc lagInteractionvars `lagInteractionvars' l_gen`var'
			}
			if "`method'" == "wrdg"{ /* WRDG */
				if !mi("`overidentification'") | !mi("`gmm'") { // in case of overidentification, or in case the user specifies a preference for GMM, estimate a system GMM model, else go with the linear IV model
					if !mi("`overidentification'"){
						loc overidentification `lagfreeVars' `lagInteractionvars'
					}
					loc eq1counter = 0
					foreach var of varlist `interactionvars'{
						loc ++eq1counter
						loc eq1vars `eq1vars' -{xb`eq1counter'}*`var'
					}
					loc eq2counter = 0
					foreach var of varlist `lagInteractionvars'{
						loc ++eq2counter
						loc eq2vars `eq2vars' -{xb`eq2counter'}*`var'
					}
					`vv' qui gmm (eq1: `depvar' `gmmfit' `eq1vars' `wrdg_contr1' - {a0}) /* y - a0 - witB - xitJ - citL
						*/ (eq2: `depvar' `gmmfit' `eq2vars' `wrdg_contr2' - {a0} -{e0}), /* y - e0 - witB - xitJ - p(cit-1L) - ... - pg(cit-1L)^g
						*/ instruments(eq1: `free' `interactionvars' `control' `translogVars' `overidentification') /* Zit1 = (1,wit,xit,c0it,cit-1)
						*/ instruments(eq2: `state' `lagfreeVars' `lagInteractionvars' `control' `translogVars') /* Zit2 = (1,xit,wit-1,cit-1)
						*/ winitial(unadjusted, independent) technique(`optimizer') conv_maxiter(`maxiter') `conv_nrtol' `init_gmm'
					/// save locals for Hanses's J and p-value
					qui estat overid
					loc hans_j: di %3.2f `r(J)'
					loc hans_p: di %3.2f `r(J_p)'
				}
				else{ /* WRDG - plain */
					tempfile wrdg
					qui save `wrdg'
						qui reg `depvar' `state' `proxyGO' `interactionvars' `free' `lagfreeVars' `control' // just to take the number of observations used in the estimation
						loc realObs = `e(N)'
						qui expand 2, gen(cons2) // double the dataset in order to stack dependent and regressors
						foreach var of local interactionvars{
							qui replace `var' = l_gen`var' if cons2 == 1 // change interactionvars to lagged values
						}
						qui ivregress gmm `depvar' `state' `proxyGO' `interactionvars' `control' (`free' = `lagfreeVars') cons2, wmatrix(unadjusted) c // run the IV regression 
						/// save locals for Hanses's J and its p-value
						loc hans_j: di %3.2f `e(J)'
						mat cV = colsof(e(V)) // find the number of instruments
						scalar cV = cV[1,1]
						loc jdf = cV - `e(rank)' + 1
						loc hans_p: di %3.2f chi2tail(`jdf', `hans_j') // this is the chi2 p-value of Hansen statistics
					qui use `wrdg', clear
					eret loc N `realObs' // post the real number of obseravtions used in estimation
				}
			}
			else if "`method'" == "mr"{ /* MrEst */
				if !mi("`overidentification'"){ // interactions are valid instruments, too
					loc overidentification  `lagInteractionvars' 
				}
				/// god forgive me: in order to overcome the difference equation in system GMM we launch both equations in level and take the initial weighting matrix
				loc instnum = (`freenum' + `statenum')
				forv i = 1/`instnum'{
					loc intVars `intVars' var_`i'
				}
				loc interinstr: list local interactionvars - intVars
				qui gmm (`depvar' `gmmfit' - {xh: `interactionvars'} `wrdg_contr1' - {a0}), quickd instruments(1: `interinstr') /*
					*/ xtinstruments(`free' `state', l(0/`lags')) winitial(xt L) onestep conv_maxiter(1)
				qui mat W1 = e(W)
				qui gmm (`depvar' `gmmfit' - {xj: `lagInteractionvars'} `wrdg_contr2' - {a0} - {e0}), quickd /*
					*/ instruments(`state' `lagInteractionvars') xtinstruments(`free' `state', l(2/`lags'))/*
					*/ /*xtinstruments(`state', l(0/`lags'))*/ winitial(xt L) onestep conv_maxiter(1)
				qui mat W2 = e(W)
				mata W_hat =  st_matrix("W1"),J(rows(st_matrix("W1")),cols(st_matrix("W2")),0) \ J(rows(st_matrix("W2")),cols(st_matrix("W1")),0), st_matrix("W2")
				mata st_matrix("W_hat", W_hat)
				/* generate the model for both equations - same parameter for c_{i,t} and c_{i,t-1} in WRDG terms */
				loc eq1counter = 0
				foreach var of varlist `interactionvars' `control'{
					loc ++eq1counter
					loc eq1vars `eq1vars' "-{xb`eq1counter'}*`var'"
				}
				loc eq2counter = 0
				foreach var of varlist `lagInteractionvars' `control'{
					loc ++eq2counter
					loc eq2vars `eq2vars' "-{xb`eq2counter'}*`var'"
				}
				/// launch the gmm with the winitial built above
				`vv' qui gmm (1: `depvar' `gmmfit' `eq1vars' /*`wrdg_contr1'*/ - {a0}) /* y - a0 - witB - xitJ - citL
					*/ (2: `depvar' `gmmfit' `eq2vars' /*`wrdg_contr2'*/ - {a0} - {e0}),/* y - e0 - witB - xitJ - p(cit-1L) - ... - pg(cit-1L)^g, where p = g = 1
					*/ instruments(1: `interinstr' `overidentification' `control') xtinstruments(1: `free' `state', l(0/`lags')) /*  
					*/ xtinstruments(2: `free' `state', l(2/`lags')) instruments(2: `state' `lagInteractionvars' `control') /* 
					*/ onestep winitial("W_hat") nocommonesample quickd technique(`optimizer') conv_maxiter(`maxiter') `conv_nrtol' `init_gmm'
				/// save locals for Hanses's J and p-value
				qui estat overid
				loc hans_j: di %3.2f `r(J)'
				loc hans_p: di %3.2f `r(J_p)'
			}
			else{ /* Robinson / ACF */ 
				qui ivregress gmm `depvar' `state' `proxyGO' `control' `lagInteractionvars' (`free' = `lagfreeVars'), vce(cluster `id') /* this is working! */
				/// save locals for Hanses's J and its p-value
				loc hans_j: di %3.2f `e(J)'
				mat cV = colsof(e(V)) // find the number of instruments
				scalar cV = cV[1,1]
				loc jdf = cV - `e(rank)' + 1
				loc hans_p: di %3.2f chi2tail(`jdf', `hans_j') // this is the chi2 p-value of Hansen statistics
			}
			/// save elements for result posting
			loc nObs = `e(N)'
			loc numInstr1: word count `e(inst_1)'
			loc numInstr2: word count `e(inst_2)'
			mat `__b' = e(b)
			mat `__b' = `__b'[1...,1..`colnum']
			mat `__V' = e(V)
			mat `__V' = `__V'[1..`colnum',1..`colnum']
			/// save locals for Hanses's J and its p-value
			loc hans_j: di %3.2f `e(J)'
			mat cV = colsof(e(V)) // find the number of instruments
			scalar cV = cV[1,1]
			loc jdf = cV - `e(rank)' + 1
			loc hans_p: di %3.2f chi2tail(`jdf', `hans_j') // this is the chi2 p-value of Hansen statistics
			continue, break
		}
		else{ /* if it's not WRDG or MrEst */
			if !mi("`attrition'"){
				foreach var in `interactionvars'{
					qui g l_gen`var' = l.`var'
					loc lagInterVars `lagInterVars' l_gen`var'
				}
				qui cap logit `exit' `lagInterVars'
				if _rc == 0{
					qui predict Pr_hat if e(sample), pr
				}
				else{
					di as error "No ID exits the sample. Running the estimation with no attrition"
					loc attrition ""
					qui g Pr_hat = .
				}
			}
			else{
				qui g Pr_hat = .
			}
			/// in case of ACF we don't want free variables to appear twice in the regression		
			if !mi("`acf'"){
				forv i = 1/`freenum'{
					loc freeVars `freeVars' var_`i'
				}
				loc interactionvars: list local interactionvars - freeVars
				/// generate the needed variables for translog: a local with all power of 2 + interactions and instruments - lag of free/proxy interacted with state
				if !mi("`translog'"){
					loc transNum: word count `transVars'
					forv i=1/`transNum'{
						forv j=`i'/`transNum'{
							loc ivar `: word `i' of `transVars''
							loc jvar `: word `j' of `transVars''
							loc tvarname "`ivar'X`jvar'"
							loc transNames `transNames' `tvarname'
							loc interactionTransVars `interactionTransVars' var_`i'_`j'
							qui g lagTrans_`i'_`j' = l.var_`i'_`j'
							loc lagInteractionTransVars `lagInteractionTransVars' lagTrans_`i'_`j'
						}
					}
					foreach fvar in `free' `proxyGO'{
						tempvar d`fvar'
						qui g `d`fvar'' = `fvar'^2
						qui g lagInstr`fvar' = l.`d`fvar''
						loc instrumentTransVars `instrumentTransVars' lagInstr`fvar'
						foreach svar in `state'{
							qui g instr_`fvar'`svar' = lagInstr`fvar'*`svar'
							loc instrumentTransVars `instrumentTransVars' instr_`fvar'`svar'
						}
					}
					foreach svar in `state'{
						tempvar d`svar'
						qui g `d`svar'' = `svar'^2
						qui g lagInstr`svar' = l.`d`svar''
						loc instrumentTransVars `instrumentTransVars' lagInstr`svar'
					}
					loc colnum: word count `free' `state' `control' `proxyGO' `interactionTransVars'
					mat `__b' = J(`reps',`colnum',.)
				}
			}
			loc regvars `free' `control' `interactionvars'
			loc firstRegNum: word count `regvars'
			loc regNum: word count `free' `state' `control' `proxyGO' `transVars'
			
			/// first stage
			qui _xt, trequired
			qui reg `depvar' `regvars'
			/// generating "freeFit" as the fit of free variables to be subtracted to Y
			tempvar freeFit
			qui g `freeFit' = 0 
			foreach var in `free'{
				scalar b_`var' = _b[`var']
				qui replace `freeFit' = `freeFit' + (b_`var'*`var')
			}
			mat `firstb' = e(b)
			mat `robV' = e(V)
			qui predict phihat if e(sample), xb
			qui g phihat_lag = l.phihat
			
			/// retrieve starting points through an OLS estimation --> first round only, not during bootstrap
			qui reg `depvar' `toLagVars' `interactionTransVars' if `touse' == 1
			if !mi("`init'"){
				mat ols_s = `init' // THIS PART IS JUST A TRYOUT IN ORDER TO MAKE THE COMMAND WORK FOR OUR PURPOSES
			}
			else{
				mat tmp = e(b)
				mata: st_matrix("ols_s",st_matrix("tmp"):+ rnormal(1,1,0,.01)) // we add some "noise" to OLS results in order not to block the optimizer
			}
			/// save the first stage results 
			if !mi("`fsresiduals'"){
				tempfile fsres
				tempvar FSfit
				mat score `FSfit' = `firstb'
				qui g `fsresiduals' = `depvar' - `FSfit'
				qui save `fsres'
			}
			qui g res = 0
			if mi("`acf'"){ /* OP and LP (non-corrected) second stage */ 
				/// here we generate a tempvar with the fitted value of all the free variables 
				qui replace phihat = phihat - `freeFit'
				qui replace res = `depvar' - `freeFit'
				qui replace phihat_lag = l.phihat
				mat init = ols_s[1...,(`freenum'+1)..`regNum']'
				loc toLagVars `state' `control' `proxyGO'
				loc laggedVars: list local laggedVars - lagfreeVars
				loc instrumentVars `state' `control' `lproxyGO'
				/// here we launch the mata routine for OP or LP 
				foreach var of varlist `toLagVars' `laggedVars' /*`lagfreeVars'*/ phihat_lag phihat{
					qui drop if mi(`var')
				}
				/// the routine cannot fail the first estimation - otherwise it would fail with the point estimates. We capture the boot repetitions
				qui mata: opt_mata(st_matrix("init"),&foplp(),"`optimizer'","phihat","phihat_lag",/*
					*/"`toLagVars'","`laggedVars'","`touse'",`maxiter',`tolerance',"`evaluator'","Pr_hat","res","`instrumentVars'","`endogenous'")
				mat `__b'[1,1] = `firstb'[1...,1..(`freenum')],r(betas)
			}
			else{ /* ACF second stage */ 
				if !mi("`overidentification'"){ // lag intereactions are valid instruments for the first equation, too
					loc overidentification `lagInteractionvars'
				}
				loc toLagVars `free' `state' `control' `proxyGO' `interactionTransVars'
				loc laggedVars `lagfreeVars' `lagstateVars' `lagcontrolVars' `lproxyGO' `lagInteractionTransVars'
				loc instrumentVars `lagfreeVars' `state' `control' `lproxyGO' `instrumentTransVars' `overidentification'
				/// here we launch the mata routine for ACF
				foreach var of varlist `laggedVars' `lagfreeVars' phihat_lag{
					qui drop if mi(`var')
				}
				loc betaNum: word count `free' `state' `proxyGO' `control' `interactionTransVars'
				mat init = ols_s[1...,1..(`betaNum')]'
				qui mata: opt_mata(st_matrix("init"),&facf(),"`optimizer'","phihat","phihat_lag",/*
					*/ "`toLagVars'","`laggedVars'","`touse'",`maxiter',`tolerance',"`evaluator'","Pr_hat","res","`instrumentVars'","`endogenous'")
				mat `__b'[1,1] = r(betas)
			}
		}
		
		if !inlist("`method'","wrdg","mr","rob"){
			clear 
			/// generate the varCovar matrix for the bootstrapped estimates
			qui svmat `__b'
			mat `__b' = `__b'[1,1...] 
		}
	
	restore 
	
	/// merge the first stage residuals
	if !mi("`fsresiduals'"){
		qui merge 1:1 `id' `t' using `fsres', nogen keepusing(`fsresiduals')
	}
	
	mat coleq `__b' = ""
	mat colnames `__b' = `free' `state' `control' `proxyGO' `transNames'
	
	/// Display results - ereturn
	eret clear
	
	eret post `__b' 
	
	if !mi("`acf'"){
		loc correction "ACF corrected"
	}
	if !mi("`fsresiduals'"){
		eret loc FSres "`fsresiduals'"
	}
	
	eret loc cmd "prodest"
	eret loc depvar "`depvar'"
	eret loc free "`free'"
	eret loc state "`state'"
	eret loc proxy "`proxy'"
	eret loc controls "`control'"
	eret loc endogenous "`endogenous'"
	eret loc method "`strMethod'"
	eret loc model "`model'"
	eret loc technique "`optimizer'"
	eret loc idvar "`id'"
	eret loc timevar "`t'"
	eret loc correction "`correction'"
	eret loc predict "prodest_p_m"
	eret loc PFtype "`PFtype'"
	eret loc gmm "`gmm'"
	
end program


//*---------------------------------------------------------------------*/

/// defining mata routines for optimization
capture mata mata drop opt_mata() 
capture mata mata drop facf()
capture mata mata drop foplp() 


*version 7
mata:

/*---------------------------------------------------------------------*/

	void foplp(todo,betas,X,lX,PHI,LPHI,RES,Z,PR_HAT,ENDO,crit,g,H)
	{
		OMEGA = PHI-X*betas'
		OMEGA_lag = LPHI-lX*betas'
		OMEGA_lag2 = OMEGA_lag:*OMEGA_lag
		OMEGA_lag3 = OMEGA_lag2:*OMEGA_lag
		/* IF clause in order to see whether we have to use the "exit" variable */
		if (!missing(PR_HAT)){
			PR_HAT2 = PR_HAT:*PR_HAT
			PR_HAT3 = PR_HAT2:*PR_HAT
			OMEGA_lag_pol = (J(rows(PHI),1,1),OMEGA_lag,OMEGA_lag2,OMEGA_lag3,PR_HAT,PR_HAT2,PR_HAT3,PR_HAT:*OMEGA_lag,PR_HAT2:*OMEGA_lag,PR_HAT:*OMEGA_lag2,ENDO)
		}
		else{
			OMEGA_lag_pol = (J(rows(PHI),1,1),OMEGA_lag,OMEGA_lag2,OMEGA_lag3,ENDO)
		}
		g_b = invsym(OMEGA_lag_pol'OMEGA_lag_pol)*OMEGA_lag_pol'OMEGA
		XI = RES-X*betas'-OMEGA_lag_pol*g_b
		crit = (XI)'*(XI)
	}
	
/*---------------------------------------------------------------------*/
	
	void facf(todo,betas,X,lX,PHI,LPHI,RES,Z,PR_HAT,ENDO,crit,g,H)
	{
		W = invsym(Z'Z)/(rows(Z))
		OMEGA = PHI-X*betas'
		OMEGA_lag = LPHI-lX*betas'
		OMEGA_lag2 = OMEGA_lag:*OMEGA_lag
		OMEGA_lag3 = OMEGA_lag2:*OMEGA_lag
		/* IF clause in order to see whether we have to use the "exit" variable */
		if (!missing(PR_HAT)){
			PR_HAT2 = PR_HAT:*PR_HAT
			PR_HAT3 = PR_HAT2:*PR_HAT
			OMEGA_lag_pol = (J(rows(PHI),1,1),OMEGA_lag,OMEGA_lag2,OMEGA_lag3,PR_HAT,PR_HAT2,PR_HAT3,PR_HAT:*OMEGA_lag,PR_HAT2:*OMEGA_lag,PR_HAT:*OMEGA_lag2,ENDO)
		}
		else{
			OMEGA_lag_pol = (J(rows(PHI),1,1),OMEGA_lag,OMEGA_lag2,OMEGA_lag3,ENDO)
		}
		g_b = invsym(OMEGA_lag_pol'OMEGA_lag_pol)*OMEGA_lag_pol'OMEGA
		XI = OMEGA-OMEGA_lag_pol*g_b
		crit = (Z'XI)'*W*(Z'XI)
	}
	
/*---------------------------------------------------------------------*/
	
void opt_mata(init, f, opt, phi, lphi, tolag, lagged, touse, maxiter, tol, eval, | Pr_hat, res, instr, endogenous)
	{
		st_view(RES=.,.,st_tsrevar(tokens(res)), touse)
		st_view(PHI=.,.,st_tsrevar(tokens(phi)), touse)
		st_view(LPHI=.,.,st_tsrevar(tokens(lphi)), touse)
		st_view(Z=.,.,st_tsrevar(tokens(instr)), touse)
		st_view(X=.,.,st_tsrevar(tokens(tolag)), touse)
		st_view(lX=.,.,st_tsrevar(tokens(lagged)), touse)
		st_view(PR_HAT=.,.,st_tsrevar(tokens(Pr_hat)), touse)
		st_view(ENDO=.,.,st_tsrevar(tokens(endogenous)), touse) 
		S = optimize_init()
		optimize_init_argument(S, 1, X)
		optimize_init_argument(S, 2, lX)
		optimize_init_argument(S, 3, PHI)
		optimize_init_argument(S, 4, LPHI)
		optimize_init_argument(S, 5, RES)
		optimize_init_argument(S, 6, Z)
		optimize_init_argument(S, 7, PR_HAT)
		optimize_init_argument(S, 8, ENDO)
		optimize_init_evaluator(S, f)
		/*optimize_init_evaluatortype(S, "d0")*/
		optimize_init_evaluatortype(S, eval)
		optimize_init_conv_maxiter(S, maxiter)
		optimize_init_conv_nrtol(S, tol)
		/*optimize_init_evaluatortype(S, "gf0")*/
		optimize_init_technique(S, opt)
		optimize_init_nmsimplexdeltas(S, 0.00001)
		optimize_init_which(S,"min")
		optimize_init_params(S,init')
		p = optimize(S)
		st_matrix("r(betas)", p)
		
		/* tryout to find Newey-type Cov estimator 
		st_matrix("r(gradient)", optimize_result_gradient(S))
		st_matrix("r(score)", optimize_result_scores(S))
		Hl = optimize_result_gradient(S)
		V = optimize_result_V(S)
		Vhh = invsym(Hl)*V*invsym(Hl)'
		st_matrix("r(Vhh)", Vhh)*/
	}
	
/*---------------------------------------------------------------------*/
	
end 



