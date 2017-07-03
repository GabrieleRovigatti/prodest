{smcl}
{* *! version 1.0.1  15sep2016}{...}
{* *! version 1.0.2  21sep2016}{...}

{cmd:help prodest}{right:also see:  {help prodest predict:prodest_predict}}
{hline}

{title:Title}

{p2colset 5 14 21 2}{...}
{p2col:{hi:prodest} {hline 1} Production Function Estimation}{p_end}
{p2colreset}{...}

{title:Syntax}

{phang}
Olley and Pakes (OP, 1996) methodology

{p 8 16 2}{cmd:prodest} {it:depvar} {ifin}
, free(varlist) proxy(varlist) state(varlist) {cmdab:met:hod(op)} [control(varlist) acf id(varname) t(varname) reps(#)
{cmdab:va:lueadded} {cmdab:l:evel(#)} poly(#) seed(#) {cmdab:fsres:idual(newname)} {cmdab:endo:genous(varlist)} {cmdab:trans:log} {it:{help prodest##optoptions:opt_options}} {it:{help prodest##opoptions:OP_options}}]

{phang}
Levinsohn and Petrin (LP, 2003) methodology

{p 8 16 2}{cmd:prodest} {it:depvar} {ifin}
, free(varlist) proxy(varlist) state(varlist) {cmdab:met:hod(lp)} [control(varlist) {cmdab:acf} id(varname) t(varname) reps(#)
{cmdab:va:lueadded} {cmdab:l:evel(#)} poly(#) seed(#) {cmdab:fsres:idual(newname)} {cmdab:endo:genous(varlist)} {cmdab:trans:log} {it:{help prodest##optoptions:opt_options}} {it:{help prodest##lpoptions:LP_options}}]

{phang}
Wooldridge (WRDG, 2009) methodology

{p 8 16 2}{cmd:prodest} {it:depvar} {ifin}
, free(varlist) proxy(varlist) state(varlist) {cmdab:met:hod(wrdg)} [control(varlist) {cmdab:acf} id(varname) t(varname) reps(#)
{cmdab:va:lueadded} {cmdab:l:evel(#)} poly(#) seed(#) {cmdab:fsres:idual(newname)} lags(#) {cmdab:endo:genous(varlist)} {it:{help prodest##optoptions:opt_options}} {it:{help prodest##wrdgoptions:WRDG_options}}]

{phang}
Mollisi and Rovigatti (MrEst, 2017) methodology

{p 8 16 2}{cmd:prodest} {it:depvar} {ifin}
, free(varlist) proxy(varlist) state(varlist) {cmdab:met:hod(mr)} [control(varlist) {cmdab:acf} id(varname) t(varname) reps(#)
{cmdab:va:lueadded} {cmdab:l:evel(#)} poly(#) seed(#) {cmdab:fsres:idual(newname)} {cmdab:endo:genous(varlist)} {it:{help prodest##optoptions:opt_options}} {it:{help prodest##mroptions:Mr_options}}]

{synoptline}
{p 4 6 2}{cmdab:met:hod(name)} accepts values {it:{help prodest##opoptions:op}} (Olley-Pakes), {it:{help prodest##lpoptions:lp}} (Levinshon-Petrin), {it:{help prodest##wrdgoptions:wrdg}} (Wooldridge) and {it:{help prodest##mroptions:mr}} (MrEst). Default is {it:{help prodest##lpoptions:lp}}{p_end}
{p 4 6 2}A panel and a time variable must be specified. Use {helpb xtset} or the {opt id(varname)} and {opt t(varname)} options.{p_end}

{synoptline}

{marker opoptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :OP_options}
{synoptline}
{syntab:Model}
{p2coldent:* {cmdab:state(}{opt varlist)}}state variable(s). Ln(capital) in OP, 1996 {p_end}
{p2coldent:* {cmdab:free(}{opt varlist)}}free variable(s). Ln(labour) in OP, 1996{p_end}
{p2coldent:* {cmdab:proxy(}{opt varlist)}}proxy variable(s). Ln(investment) in OP, 1996{p_end}
{synopt :{cmdab:control(}{opt varlist)}}control variable(s) to be included{p_end}
{synopt :{cmdab:endo:genous(}{opt varlist)}}endogenous variable(s) to be included{p_end}
{synopt :{opt acf}}apply the Ackerberg, Caves and Frazer (2015) correction{p_end}
{synopt :{cmdab:va:lueadded}}{it:depvar} is the value added of output. Default is gross output{p_end}
{synopt :{cmdab:att:rition}}correct for attrition - i.e. firm exit - in the data{p_end}
{synopt :{cmdab:trans:log}}use a translog production function for estimation - ACF only{p_end}

{syntab:Optimization}
{synopt :{cmdab:opt:imizer(}{it:{help prodest##optimizers:opttype}})}available optimizers are Nelder-Mead (nm), modified Newton-Raphson (nr), 
Davidon-Fletcher-Powell (dfp), Broyden-Fletcher-Goldfarb-Shanno (bfgs) and Berndt-Hall-Hall-Hausman (bhhh){p_end}
{synopt :{cmdab:max:iter}}max number of iterations, default is 10,000{p_end}

{syntab:Other}
{synopt :{cmdab:id(}{opt varname})}{it:panelvar} to {helpb xtset} the data{p_end}
{synopt :{cmdab:t(}{opt varname})}{it:timevar} to {helpb xtset} the data{p_end}
{synopt :{opt reps(#)}}number of bootstrap repetitions; minimum is 2, default is 5{p_end}
{synopt :{opt poly(#)}}degree of polynomial approximations; ranges 3 to 6, default is 3{p_end}
{synopt :{opt seed(#)}}seed to be set; integer, default is 12345{p_end}
{synopt :{cmdab:fsres:iduals(}{opt newname})}save first stage residuals in var {it:newname}{p_end}

{syntab:Reporting}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}.{p_end}

{synoptline}

{marker lpoptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :LP_options}
{synoptline}
{syntab:Model}
{p2coldent:* {cmdab:state(}{opt varlist)}}state variable(s). Ln(capital) in LP, 2003 {p_end}
{p2coldent:* {cmdab:free(}{opt varlist)}}free variable(s). Ln(labour) in LP, 2003{p_end}
{p2coldent:* {cmdab:proxy(}{opt varlist)}}proxy variable(s). Ln(intermediate inputs) in LP, 2003{p_end}
{synopt :{cmdab:control(}{opt varlist)}}control variable(s) to be included{p_end}
{synopt :{cmdab:endo:genous(}{opt varlist)}}endogenous variable(s) to be included{p_end}
{synopt :{opt acf}}apply the Ackerberg, Caves and Frazer (2015) correction{p_end}
{synopt :{cmdab:va:lueadded}}{it:depvar} is the value added of output. Default is gross output{p_end}
{synopt :{cmdab:att:rition}}correct for attrition - i.e. firm exit - in the data{p_end}
{synopt :{cmdab:trans:log}}use a translog production function for estimation - ACF only{p_end}

{syntab:Optimization}
{synopt :{cmdab:opt:imizer(}{it:{help prodest##optimizers:opttype}})}available optimizers are Nelder-Mead (nm), modified Newton-Raphson (nr), 
Davidon-Fletcher-Powell (dfp), Broyden-Fletcher-Goldfarb-Shanno (bfgs) and Berndt-Hall-Hall-Hausman (bhhh){p_end}
{synopt :{cmdab:max:iter}}max number of iterations, default is 10,000{p_end}

{syntab:Other}
{synopt :{cmdab:id(}{opt varname})}{it:panelvar} to {helpb xtset} the data{p_end}
{synopt :{cmdab:t(}{opt varname})}{it:timevar} to {helpb xtset} the data{p_end}
{synopt :{opt reps(#)}}number of bootstrap repetitions; minimum is 2, default is 5{p_end}
{synopt :{opt poly(#)}}degree of polynomial approximations; ranges 3 to 6, default is 3{p_end}
{synopt :{opt seed(#)}}seed to be set; integer, default is 12345{p_end}
{synopt :{cmdab:fsres:iduals(}{opt newname})}save first stage residuals in var {it:newname}{p_end}

{syntab:Reporting}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}.{p_end}

{synoptline}

{marker wrdgoptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :WRDG_options}
{synoptline}
{syntab:Model}
{p2coldent:* {cmdab:state(}{opt varlist)}}state variable(s). Ln(capital) in WRDG, 2009{p_end}
{p2coldent:* {cmdab:free(}{opt varlist)}}free variable(s). Ln(labour) in WRDG, 2009{p_end}
{p2coldent:* {cmdab:proxy(}{opt varlist)}}proxy variable(s). Ln(intermediate inputs) in WRDG, 2009{p_end}
{synopt :{cmdab:control(}{opt varlist)}}control variable(s) to be included{p_end}
{synopt :{cmdab:endo:genous(}{opt varlist)}}endogenous variable(s) to be included{p_end}
{synopt :{cmdab:va:lueadded}}{it:depvar} is the value added of output. Default is gross output{p_end}

{syntab:Optimization}
{synopt :{cmdab:opt:imizer(}{it:{help prodest##optimizers:opttype}})}available optimizers are Nelder-Mead (nm), modified Newton-Raphson (nr), 
Davidon-Fletcher-Powell (dfp), Broyden-Fletcher-Goldfarb-Shanno (bfgs) and Berndt-Hall-Hall-Hausman (bhhh){p_end}
{synopt :{cmdab:max:iter}}max number of iterations, default is 10,000{p_end}

{syntab:Other}
{synopt :{cmdab:id(}{opt varname})}{it:panelvar} to {helpb xtset} the data{p_end}
{synopt :{cmdab:t(}{opt varname})}{it:timevar} to {helpb xtset} the data{p_end}
{synopt :{opt poly(#)}}degree of polynomial approximations; ranges 3 to 6, default is 3{p_end}
{synopt :{opt seed(#)}}seed to be set; integer, default is 12345{p_end}
{synopt :{cmdab:fsres:iduals(}{opt newname})}save first stage residuals in var {it:newname}{p_end}

{syntab:Reporting}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}.{p_end}

{synoptline}

{marker mroptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :Mr_options}
{synoptline}
{syntab:Model}
{p2coldent:* {cmdab:state(}{opt varlist)}}state variable(s){p_end}
{p2coldent:* {cmdab:free(}{opt varlist)}}free variable(s){p_end}
{p2coldent:* {cmdab:proxy(}{opt varlist)}}proxy variable(s){p_end}
{synopt :{cmdab:control(}{opt varlist)}}control variable(s) to be included{p_end}
{synopt :{cmdab:endo:genous(}{opt varlist)}}endogenous variable(s) to be included{p_end}
{synopt :{cmdab:va:lueadded}}{it:depvar} is the value added of output. Default is gross output{p_end}
{synopt :{opt lags(#)}}lags to be used as Blundell-Bond-type instruments{p_end}

{syntab:Optimization}
{synopt :{cmdab:opt:imizer(}{it:{help prodest##optimizers:opttype}})}available optimizers are Nelder-Mead (nm), modified Newton-Raphson (nr), 
Davidon-Fletcher-Powell (dfp), Broyden-Fletcher-Goldfarb-Shanno (bfgs) and Berndt-Hall-Hall-Hausman (bhhh){p_end}
{synopt :{cmdab:max:iter}}max number of iterations, default is 10,000{p_end}

{syntab:Other}
{synopt :{cmdab:id(}{opt varname})}{it:panelvar} to {helpb xtset} the data{p_end}
{synopt :{cmdab:t(}{opt varname})}{it:timevar} to {helpb xtset} the data{p_end}
{synopt :{opt poly(#)}}degree of polynomial approximations; ranges 3 to 6, default is 3{p_end}
{synopt :{opt seed(#)}}seed to be set; integer, default is 12345{p_end}
{synopt :{cmdab:fsres:iduals(}{opt newname})}save first stage residuals in var {it:newname}{p_end}

{syntab:Reporting}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}.{p_end}

{synoptline}

{marker optoptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :Optimization Options}
{synoptline}
{syntab:Optimization}
{synopt :{cmdab:opt:imizer}({it:{help prodest##optimizers:opttype}})}available optimizers are Gauss-Newton (gn), modified Newton-Raphson (nr), 
Davidon-Fletcher-Powell (dfp) and Broyden-Fletcher-Goldfarb-Shanno (bfgs){p_end}
{synopt :{cmdab:max:iter}(#)}max number of iterations, default is 10,000{p_end}
{synopt :{cmdab:eval:uator}({opt string})}{it:{help mf_optimize##i_evaluator:init_evaluator_type()}}, default is d0 (gf0 for bhhh){p_end}
{synopt :{cmdab:tol:erance(#)}}{it:{help mf_optimize##i_ptol:init_conv_nrtol()}}, default is e-05{p_end}


{marker optimizers}{...}
{synoptset 35 tabbed}{...}
{synopthdr :Optimizer Techniques}
{synoptline}
{synopt :{opt nm}}Nelder-Mead*{p_end}
{synopt :{opt nr}}modified Newton-Raphson{p_end}
{synopt :{opt dfp}}Davidon-Fletcher-Powell{p_end}
{synopt :{opt bfgs}}Broyden-Fletcher-Goldfarb-Shanno{p_end}
{synopt :{opt gn}}Gauss-Newton**{p_end}
{synopt :{opt bhhh}}Berndt-Hall-Hall-Hausman*{p_end}
{synoptline}
	*op and lp methods only
	**wrdg and mr methods only
	
{title: Description}

{pstd}
{cmd:prodest} aims at estimating production functions using the control function approach. It includes Olley-Pakes (OP), Levinshon-Petrin (LP), Wooldridge (WRDG) and Ackerberg-Caves-Frazer (ACF) estimation methodologies. A new technique (MrEst) has been added in order to deal with short panels.
By default, the command requires the log gross output [or value added] variable (y_i,t), a set of free variables (typically log labour, w_i,t), a set of state variables (typically log capital, k_i,t) and a set of proxy variables. 

{pstd}
Consider the following Cobb-Douglas production technology for firm i at time t:

		y_i,t = alpha + w_i,t*beta + k_i,t*gamma + omega_i,t + epsilon_i,t

{pstd}
where y_it is the (log) gross output, w_it is a 1xJ vector of (log) free variables, k_it is a 1xK vector of state variables and epsilon_it is a normally distributed idiosyncratic error term. 
The random component omega_it is the unobserved technical efficiency parameter. It evolves according to a first-order Markov process:

		omega_it = E(omega_i,t | omega_i,t-1) + u_i,t = g(omega_i,t-1) + u_i,t

{pstd}	
and u_i,t is a random shock component assumed to be uncorrelated with the technical 
efficiency, the state variables in k_it and the lagged free variables (w_i,t-1). 
Productivity estimation in OP and LP methods (and their ACF corrections) is performed in two steps: 

{cmd:i)} {pstd} the OP method relies on the following set of assumptions:

	{cmd:a)} inv_i,t = inv(k_i,t , omega_i,t) - investments are a function of both the state variable and the technical efficiency parameter;
	{cmd:b)} inv_i,t is strictly monotone in omega_i,t ;
	{cmd:c)} omega_i,t is scalar unobservable in i_i,t = i(.) ;
	{cmd:d)} the levels of inv_i,t and k_i,t are decided at time t-1; the level of the free variable, w_i,t , is decided after the shock u_i,t realizes. 

{pstd}
Assumptions a)-d) ensure the invertibility of inv_i,t in omega_i,t and lead to the partially identified model:

	y_i,t = alpha + w_i,t*beta + k_i,t*gamma + h(inv_i,t, k_i,t ) + epsilon_i,t = alpha + w_i,t*beta + psi(inv_i,t , k_i,t) + epsilon_i,t
	
{pstd}
which can be estimated by non parametric approach - First Stage. Exploiting the Markovian nature of the productivity process one can exploit assumption d) as moment conditions to estimate the production function parameters - Second stage.
Exploting the resisual e_i,t of:

	y_i,t - w_i,t*beta_hat = alpha + k_i,t*gamma + g(omega_i,t-1,chi_i,t) + e_i,t
	
{pstd}
where g(.) is typically left unspecified and approximated by a n-th order polynomial and chi_i,t is an inidicator function for the attrition in the market.

{cmd:ii)} 
{pstd} LP aimed to overcome the empirical issue of zeros in the investment data: They proposed instead to use intermediate inputs as a proxy 
variable for omega_i,t , under the following set of assumptions:

	{cmd:a)} firms immediately adjust the level of inputs according to demand function m(omega_i,t , k_i,t) after the technical efficiency shock is realized;
	{cmd:b)} m_i,t is strictly monotone in omega_i,t ;
	{cmd:c)} omega_i,t is scalar unobservable in m_i,t = m(.) ;
	{cmd:d)} the levels of k_i,t are decided at time t-1; the level of the free variable, w_i,t , is decided after the shock u_i,t is realized.
	
{pstd}
Assumptions a)-d) ensure the invertibility of m_i,t in omega_i,t and lead to the partially identified model:
	
	y_i,t = alpha + w_i,t*beta + psi(m_i,t , k_i,t) + v_i,t

{pstd}
which can be estimated by non parametric approach - First Stage. Exploiting the Markovian nature of the productivity process one can exploit assumption d) as moment conditions to estimate the production function parameters - Second stage.
Exploting the resisual v_i,t of:

	y_i,t - w_i,t*beta_hat = alpha + k_i,t*gamma + g(omega_i,t-1,chi_i,t) + v_i,t
	
{pstd}
where g(.) is typically left unspecified and approximated by a n-th order polynomial and chi_i,t is an inidicator function for the attrition in the market.


{cmd:iii)}
{pstd} Labour demand and the control function are partially collinear. ACF proposed an alternative estimation algorithm based on the following assumptions:

	{cmd:a)} p_i,t = p(k_i,t , l_i,t , omega_i,t) is the proxy variable policy function ;
	{cmd:b)} Strict monotonicity holds for p_i,t relative to omega_i,t ;
	{cmd:c)} omega_i,t is scalar unobservable in p_i,t = p(.) ;
	{cmd:d)} The state variable are decided at time t-1. The less variable labor input, l_i,t , is chosen at t-b, where 0 < b < 1. 
		The free variables, w_i,t, are chosen in t when the firm productivity shock is realized.

{pstd}
Under this set of assumptions, the first stage is meant to remove the shock epsilon_i,t from the the output, y_i,t. As before, the policy function can replace, once inverted, the productivity term, omega_i,t, in the production function, yielding:

	y_i,t = k_i,t*gamma + w_i,t*beta + mu*l_i,t + h(p_i,t , k_i,t ,w_i,t , l_i,t) + epsilon_i,t

{pstd}
which can be estimated by non parametric approach - First Stage. Exploiting the Markovian nature of the productivity process one can exploit assumption d) as moment conditions to estimate the production function parameters - Second stage.
Exploting the resisual u_i,t of:

	omega_i,t = E(omega_i,t | omega_i,t-1) + u_i,t} = g(omega_i,t-1) + u_i,t
 

{cmd:iv)}
{pstd} Wooldridge proposed a more efficient approach to implement OP and LP methodology. The two stages are jointly estimated and the estimator does not suffer from the collinearity issue identified in OP/LP methodology.
Wooldridge's system gmm is based on the following assumptions:

	{cmd:a)} omega_i,t = g(x_i,t , p_i,t), productivity is an unknown function g(.) of state and a vector of proxy variables, p_i,t ;
	{cmd:b)} E(omega_i,t | omega_i,t-1 ,...,,omega_i,t-T) = E(omega_i,t|omega_i,t-1), t = 2,,3,...,,T. Assumption b restrains the productivity's 
		dynamic to a first order Markov chain process. 
	{cmd:c)} E(omega_i,t | omega_i,t-1)=f[omega_i,t-1], productivity is an unknown function f[.] of lagged productivity, omega_i,t-1.
 
{pstd}Under the above set of assumptions, It is possible to construct a system gmm using the  vector of residuals from
             
			y_i,t - alpha - w_i,t*beta - x_i,t*gamma - g(x_i,t , p_i,t)
	r_i,t = 
			y_i,t - alpha - w_i,t*beta - x_i,t*gamma - f[g(x_i,t-1 , p_i,t-1)]
			
{pstd}where the unknown function f(.) is approximeted by a n-th polynomial and  g(x_i,t , m_i,t) = lambda_0 + c(x_i,t , m_i,t)*lambda. In particular, g(x_i,t , m_i,t) is a linear combination of functions in (x_i,t , m_i,t)
where c_i,t are the addends of this linear combination. The residuals r_i,t are used to set the moment conditions 

	E(Z_i,t*r_i,t) =0

{pstd}with the following set of instruments:

		    z1_i,t = (1 , w_i,t , x_i,t , c_i,t)
	Z_i,t =	
		    z2_i,t = (w_i,t-1 , x_i,t , c_i,t)
	
{cmd:v)}
{pstd}Previous lags of free and state variables are potentially valid instruments in a Wooldridge-type estimation framework. However, adding instruments in such context would lead to a reduced sample size, and this could be problematic given the "large N, small T" nature of most dataset used in the related literature.
Introducing dynamic panel instruments à la Blundell-Bond is a solution to the issue: this allows to exploit the additional information in the lagged instruments without losing observations and estimation power. Defining the residual function matrix as:
			 															
			y_i,2 - alpha - w_i,2*beta - x_i,2*gamma - g(x_i,2 , p_i,2)		
			y_i,2 - alpha - w_i,2*beta - x_i,1*gamma - f[g(x_i,1 , p_i,1)]	
	r_i,t =					...										
						...										
			y_i,t - alpha - w_i,t*beta - x_i,t*gamma - g(x_i,t , p_i,t)		
			y_i,t - alpha - w_i,t*beta - x_i,t*gamma - f[g(x_i,t-1 , p_i,t-1)]

{pstd}for each panel i we define t - b the last available lag (i.e. b = 1 at t = 2, b = 2 at t = 3, etc.). Then, define Z_i the dynamic panel instrument matrix for each panel (we suppress the subscript i):

			z'_2 z'_3 ... z'_T 0    0    0   0   
			0    0    ... 0    S'_3 0    0   0   
	Z = 		0    0    ... 0    0    S'_4 0   0   
			0    0    ... 0    0    0   ...  S'_T
			0    0    ... 0    1    1   ...  1   
			
{pstd}where S_t is a 1xb vector consisting of [z_t-1, ... , z_t-b]. Usual moment conditions E[Z_i,t r_i,t] = 0 define the MrEst. 
	
{marker examples}{...}
{title:Examples}

{tab}{cmd:. insheet using "https://raw.githubusercontent.com/GabBrock/prodest/master/prodest.csv", names clear}}
{phang2}{it:({stata "insheet using https://raw.githubusercontent.com/GabBrock/prodest/master/prodest.csv, names clear":load data})}

{tab}{cmd:. xtset id year, y /* not run */}

{pstd}Olley and Pakes method{p_end}
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_investment) va met(op) poly(4) reps(40) id(id) t(year)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_investment) va met(op) poly(4) reps(40) id(id) t(year)":click to run})}

{tab}Ackerberg, Caves and Frazer correction
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_investment) va met(op) acf opt(nm) reps(50) id(id) t(year) fsresiduals(fs_acf_op)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_investment) va met(op) acf opt(nm) reps(50) id(id) t(year) fsres(fs_acf_op)":click to run})}

{pstd}Levinsohn and Petrin method{p_end}
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) opt(dfp) reps(50) id(id) t(year) fsresiduals(fs_lp)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) opt(dfp) reps(50) id(id) t(year) fsres(sf_lp)":click to run})}

{tab}Ackerberg, Caves and Frazer correction
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) acf reps(50) id(id) t(year)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) acf reps(50) id(id) t(year)":click to run})}

{pstd}Wooldridge method{p_end}
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(wrdg) poly(2) id(id) t(year)} 
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(wrdg) poly(2) id(id) t(year)":click to run})}

{pstd}MrEst method{p_end}
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(mr) lags(1) poly(2) id(id) t(year)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(mr) lags(1) poly(2) id(id) t(year)":click to run})}

{title:References}

{phang}
Ackerberg D. A., Caves K., and Frazer G., 2015 
Identification properties of recent production function estimators.
Econometrica, 83, (6), pp. 2411-2451.

{phang}
Levinsohn, J. and Petrin, A., 2003.
Estimating Produciton Functions Using Inputs to Control for Unobservables.
Review of Economic Studies, 70, pp. 317-341.

{phang}
Mollisi, V. and Rovigatti, G., 2017.
Theory and Practice of TFP Estimation: the Control Function Approach Using Stata.
CEIS Working Paper Series, No. 399.
{it:({stata "!cmd /c start https://papers.ssrn.com/sol3/papers2.cfm?abstract_id=2916753":click to download - Windows})}
{it:({stata "!xdg-open https://papers.ssrn.com/sol3/papers2.cfm?abstract_id=2916753":click to download - Unix})}
{it:({stata "!open https://papers.ssrn.com/sol3/papers2.cfm?abstract_id=2916753":click to download - MacOSX})}

{phang}
Olley, S. O. and Pakes, A., 1996.
The dynamics of productivity in the telecommunications equipment industry.
Econometrica, 64, pp. 1263-1297. 

{phang}
Wooldridge, J., 2009. 
On estimating firm-level production functions using proxy variables to control for unobservables.
Economics Letters, 104, pp. 112-114.

{title:Saved results}

{pstd}
{cmd:prodest} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of ID{p_end}
{synopt:{cmd:e(tmin)}}min number of periods{p_end}
{synopt:{cmd:e(tmean)}}average number of periods{p_end}
{synopt:{cmd:e(tmax)}}max number of periods{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:prodest}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(free)}}varname of the free variable(s){p_end}
{synopt:{cmd:e(state)}}varname of the state variable(s){p_end}
{synopt:{cmd:e(proxy)}}varname of the proxy variable(s){p_end}
{synopt:{cmd:e(control)}}varname of the control variable(s){p_end}
{synopt:{cmd:e(endogenous)}}varname of the endogenous variable(s){p_end}
{synopt:{cmd:e(technique)}}optimization technique{p_end}
{synopt:{cmd:e(idvar)}}ID variable{p_end}
{synopt:{cmd:e(timevar)}}time variable{p_end}
{synopt:{cmd:e(method)}}estimation methodology{p_end}
{synopt:{cmd:e(model)}}value added or gross output model{p_end}
{synopt:{cmd:e(correction)}}correction - ACF{p_end}
{synopt:{cmd:e(hans_j)}}Hansen's J - wrdg{p_end}
{synopt:{cmd:e(hans_p)}}Hansen's J p-value{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{title:Authors}

{pstd}Gabriele Rovigatti{p_end}
{pstd}University of Rome Tor Vergata{p_end}
{pstd}Einaudi Institute for Economics and Finance{p_end}
{pstd}Rome, Italy{p_end}
{pstd}gabriele.rovigatti@gmail.com{p_end}

{pstd}Vincenzo Mollisi{p_end}
{pstd}University of Rome Tor Vergata{p_end}
{pstd}Rome, Italy{p_end}
{pstd}vincenzo.mollisi@gmail.com{p_end}

{title:Also see}

{psee}
Online: {helpb opreg}, {helpb levpet}, {helpb acfest}{p_end}

{phang}
Yasar, M., Raciborski, R. and Poi, B., 2008.
Production function estimation in Stata using the Olley and Pakes method.
The Stata Journal, 8, pp. 221-231.

{phang}
Petrin, A., Poi, B. and Levinsohn, J., 2004.
Production function estimation in Stata using inputs to control for unobservables.
The Stata Journal, 4, pp. 113-123.

{phang}
Manjon, M. and Manez, J., 2016.
Production function estimation in Stata using the Ackerberg-Caves-Frazer method.
The Stata Journal, 16, pp. 900-916.

{title:Bug Reporting}

{psee}
{opt prodest} is part of an ongoing project. Thus, it might contain errors and malfunctions. Please submit bugs, comments and suggestions via email to:

{pstd}	gabriele.rovigatti@gmail.com

{psee}
Please try to include steps to reproduce issues and possibly the output of the Stata command -creturn list-.


