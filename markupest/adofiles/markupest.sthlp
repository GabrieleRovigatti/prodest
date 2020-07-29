{smcl}
{* *! version 1.0.1  25may2020}{...}

{cmd:help markupest}
{hline}

{title:Title}

{p2colset 5 14 21 2}{...}
{p2col:{hi:markupest} {hline 1} Markup Estimation}{p_end}
{p2colreset}{...}

{title:Syntax}

{syntab:MICRO}

{phang}
De Loecker and Warzynski (DLW, 2012) 

{p 8 16 2}{cmd:markupest} {it:markupvar} {ifin}
, {cmdab:input:var(varname)} {cmdab:met:hod(dlw)} [{cmdab:verb:ose} {cmdab:saveb:eta} replace {cmdab:corr:ected} id(varname) t(varname) output(varname) free(varlist) proxy(varlist) state(varname) {cmdab:prodestopt:ions}({it:{help prodest:see prodest}}) {cmdab:est:imates}({it:{help estimates_store:estimates name}})}]

{syntab:MACRO}

{phang}
Hall (2018)

{p 8 16 2}{cmd:markupest} {it:markupvar} {ifin}
, {cmdab:met:hod(hall)} [{cmdab:verb:ose} {cmdab:saveb:eta} W(varname) replace id(varname) t(varname) GO(varname) deltago(varname) pgo(varname) DELTAvars(varlist) PRICEvars(varlist) INSTRuments(varlist) timevarying]

{phang}
Roeger (1995)

{p 8 16 2}{cmd:markupest} {it:markupvar} {ifin}
, {cmdab:met:hod(roeger)} [{cmdab:verb:ose} {cmdab:saveb:eta} W(varname) replace id(varname) t(varname) GO(varname) INputs(varlist)]

{synoptline}
{p 4 6 2}{cmdab:met:hod(name)} accepts values {it:{help markupest##dlwoptions:dlw}} (De Loecker-Warzynski), {it:{help markupest##roegeroptions:roeger}} (Roeger), {it:{help markupest##halloptions:hall}} (Hall, 1988, 2018){p_end}

{synoptline}

{marker dlwoptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :DLW_options}
{synoptline}
{syntab:General}
{synopt :}A panel and a time variable must be specified. Use {helpb xtset} or the {opt id(varname)} and {opt t(varname)} options{p_end}
{p2coldent:* {cmdab:met:hod(dlw)}}estimation method{p_end}
{p2coldent:* {cmdab:inputvar(}{opt varlist)}}variable input used for markup estimation{p_end}
{synopt :{cmdab:corr:ected}}correct markups with the first-stage resdiuals. Requires PF estimation with the {helpb prodest} option {opt fsreadisuals(newvar)}{p_end}
{synopt :{cmdab:verb:ose}}print on-screen the estimates of beta{p_end}
{synopt :{cmdab:saveb:eta}}store the estimated beta in variables _b + varname {p_end}
{synopt :{cmdab:replace}}replace {opt markupvar} {p_end}

{syntab:PF estimation}
{p2coldent:+ {cmdab:output(}{opt varname)}}log() of output variable - e.g., value added - used for PF estimation{p_end}
{p2coldent:+ {cmdab:state(}{opt varname)}}log() of state variable - capital - used for PF estimation{p_end}
{p2coldent:+ {cmdab:free(}{opt varlist)}}log() of free variable(s) - e.g., labor - used for PF estimation{p_end}
{p2coldent:+ {cmdab:proxy(}{opt varlist)}}log() of proxy variable(s) - e.g. intermediate input - used for PF estimation{p_end}
{synopt :{cmdab:prodestopt:ions(string)}}all available options for PF estimation in prodest. see {helpb prodest}{p_end}

{syntab:Saved results}
{p2coldent:° {cmdab:est:imates}(name)}{helpb prodest} saved results to be used for markup estimation (see {helpb estimates_store}){p_end}

{synoptline}

{marker halloptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :HALL_options}
{synoptline}
{syntab:General}
{synopt :}A panel and a time variable must be specified. Use {helpb xtset} or the {opt id(varname)} and {opt t(varname)} options{p_end}
{p2coldent:* {cmdab:met:hod(hall)}}estimation method{p_end}
{p2coldent:* {cmdab:GO}({opt varname})}gross output{p_end}
{p2coldent:* {cmdab:deltago}({opt varname})}change in gross output{p_end}
{p2coldent:* {cmdab:pgo}({opt varname})}change in gross output prices{p_end}
{p2coldent:* {cmdab:delta:vars}({opt varlist})}changes in production inputs{p_end}
{p2coldent:* {cmdab:price:vars}({opt varlist})}changes in production input costs{p_end}
{p2coldent:* {cmdab:instr:uments}({opt varlist})}instrumental variable(s){p_end}
{synopt :{cmdab:timev:arying}}returns the time-varying version of Hall's estimator{p_end}
{synopt :{cmdab:W}({opt varname})}regression weights{p_end}
{synopt :{cmdab:verb:ose}}print on-screen the estimates of beta{p_end}
{synopt :{cmdab:saveb:eta}}store the estimated beta in variables _b + varname {p_end}

{synoptline}

{marker roegeroptions}{...}
{synoptset 29 tabbed}{...}
{synopthdr :ROEGER_options}
{synoptline}
{syntab:General}
{synopt :}A panel and a time variable must be specified. Use {helpb xtset} or the {opt id(varname)} and {opt t(varname)} options{p_end}
{p2coldent:* {cmdab:met:hod(roeger)}}estimation method{p_end}
{p2coldent:* {cmdab:go}({opt varname})}gross output{p_end}
{p2coldent:* {cmdab:in:puts}({opt varname})}production inputs - capital must be listed last{p_end}
{synopt :{cmdab:W}({opt varname})}regression weights{p_end}
{synopt :{cmdab:verb:ose}}print on-screen the estimates of beta{p_end}
{synopt :{cmdab:saveb:eta}}store the estimated beta in variables _b + varname {p_end}

{synoptline}
	
{title: Description}

{pstd}
{cmd:markupest} implements the markup estimation micro-approach proposed by De Loecker and Warzynski (AER, 2012), and the macro estimation routines proposed by Hall (JPE, 1988), Roeger (JPE, 1995) and Hall (WP, 2018).

{pstd}
{cmd:Micro.1)} De Loecker and Warzynski (2012)

{pstd}
Consider N heterogeneous firms competing à la Cournot and producing output according to a multiple-input Production Function (PF) of the form

		Y_it = F_it(I^1_it, ... , I^G_it, K_it, u_it)								(1)

where I^1_it, ... , I^G_it are variable inputs - including labor, electricity, intermediate inputs, etc. - K_it is capital and u_it is the idiosyncratic, Hicks-neutral productivity parameter.

Without loss of generality, consider a single variable input production technology in labor L_it, with an associated wage of w_it. DLW consider the dual of the profit maximization problem of the firm (i.e., the cost minimization) and setup a Lagrangian whose FOC with respect to labor reads

		dLagr^DLW / dL_it = w_it - lambda_it * dF(L_it, K_it, u_it)/dL_it 						(2)
		
note that dLagr^DLW / dY_it = lambda_it stands for the marginal cost of production at any level of output Y_it. Accordingly, the most straightforward definition of markup is the price-marginal cost ratio, i.e., mu_it = P_it / lambda_it.

Rearranging (2) and multiplying both sides by L_it / Y_it yields

		dF(L_it, K_it, u_it)/dL_it * L_it/Y_it = 1/lambda_it * w_it L_it/Y_it				(3)
		
which, in words, means that the optimal demand for labor is obtained by equalizing the output elasticity (theta^L_it = dF(L_it, K_it, u_it)/dL_it * L_it/Y_it) and an inverse function of marginal cost lambda_it and labor share (w_it L_it/Y_it).

Building on the above, the estimable formula for markup estimation reads mu_it = theta^L_it / alpha^L_it, where  alpha^L_it is the labor share on total output (i.e., input share for multiple variable input PFs). 


{pstd}
{cmd:Macro.1)} Hall (1986, 1988, 2018)

Under the assumption of Constant Returns to Scale (CSR), Hicks-neutral technology and perfect competition, the technological changes Delta log(Theta)_t = theta_t are captured by the so-called Solow residuals. Indeed, it is possible to restate output growth in terms of weighted changes in inputs (X_it) plus changes in technology:

		Delta Y_t = \sum_i alpha_it Delta X_it + theta_t								(4)
		
however, Hall (1988) showed that, when the assumption of prerfect competition does not hold, the Solow residuals do not nail down theta_t anymore, being instead a weighted sum of technological changes and growth rate of output-capital ratio (whose weights are a function of the markups. More specifically, the formula for the quantity-based Solow residual reads

		SolowResQ_it = Delta Y_it - alpha^L_t Delta L_t - alpha^M_t Delta M_t - (1 - alpha^L_t - alpha^M_t) Delta K_t = (1 - 1/mu_t) (Delta Y_t - Delta K_t) + (1/mu_t) theta_t			(5)

where (5) must be corrected for endogeneity in order to ensure a correct identificaiton. Hall (1988, 2018) proposes to use exogenous demand shifters like unexpected shocks in military spending, or changes in energy costs, as instruments. 

{pstd}
{cmd:Macro.2)} Roeger (1995)

In order to overcome the main drawbacks of the Hall model - which are i) the validity of the instruments, for all sectors, countries, or periods of time, and ii) the availability of the data. In order to do that, Roeger proposes to exploit the dual problem of (4), that is, use the cost minimization of the firm, too, and obtain a price-based Solow residual to be constrasted with its quantity-based counterpart. The former is defined as

		SolowResP_it = Delta p_t - alpha^L_t Delta w_t - alpha^M_t Delta m_t - (1 - alpha^L_t - alpha^M_t) Delta r_t = (1 - 1/mu_t) (Delta p_t - Delta r_t) - (1/mu_t) theta_t			(6)

where p_t is gross output price, w_t stands for the wage level, whereas m_t and r_t capture input and capital prices, respectively. The main intuition behind the model is that, by summing up (5) and (6), the technological changes are net out of the equation, and all resulting variables are expressed in nominal, and observable, terms. More specifically, the resulting equation is

		SolowResQ_t - SolowResP_t = (Delta p_t + Delta Y_t) - alpha^L_t(Delta w_t + Delta L_t) - alpha^M_t(Delta m_t + Delta M_t) - (1 - alpha^L_t - alpha^M_t)(Delta r_t + Delta K_t) =
		= (1 - 1/mu_t)[(Delta p_t + Delta Q_t) - (Delta r_t + Delta K_t)]
		
which is estimable with an OLS model. More specifically, if we define Y_t = SolowResQ_t - SolowResP_it, X_t = [(Delta p_t + Delta Q_t) - (Delta r_t + Delta K_t)], then Y_t = beta X_t + epsilon_t is identified only if one assumes a constant markup (i.e., mu_t = mu) which can be retrieved by inverting the estimated beta - i.e., mu^R = 1 / (1 - beta^R).
	
{marker examples}{...}
{title:Examples}

{pstd}
{cmd:Ex.1)} DLW method: markup estimation on simulated data

{tab}{cmd:. use "https://raw.githubusercontent.com/GabrieleRovigatti/prodest/master/markupest/data/f_DGP3_replica.dta", clear}
{phang2}{it:({stata "use https://raw.githubusercontent.com/GabrieleRovigatti/prodest/master/markupest/data/f_DGP3_replica.dta, clear":load data})}

{pstd}Markup estimation: Translog ACF, corrected for first stage{p_end}
{tab}{cmd:. markupest mkupD3_m1_t_corr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm1) prodestopt("poly(3) acf trans va") corrected verbose}
{phang2}{it:({stata `"markupest mkupD3_m1_t_corr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm1) prodestopt("poly(3) acf trans va") corrected verbose"':generate markups})}

{pstd}Markup estimation: Cobb-Douglas, not corrected{p_end}
{tab}{cmd:. markupest mkupD3_m1_t_corr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm1) prodestopt("poly(3) acf trans va") corrected verbose}
{phang2}{it:({stata `"markupest mkupD3_m4_cb_uncorr, method(dlw) output(lny) inputvar(lnl) free(lnl) state(lnk) proxy(lnm1) prodestopt("poly(2) acf va") verbose"':generate markups})}

{pstd}
{cmd:Ex.2)} Hall and Roeger methods: markup estimation on BLS data

{tab}From "https://www.bls.gov/mfp/special_requests/klemscombinedbymeasure.xlsx"
{tab}{cmd:. use "https://raw.githubusercontent.com/GabrieleRovigatti/prodest/master/markupest/data/bls_data.dta", clear}
{phang2}{it:({stata "use https://raw.githubusercontent.com/GabrieleRovigatti/prodest/master/markupest/data/bls_data.dta, clear":load data})}

{pstd}Hall - time invariant{p_end}
{tab}{cmd:. bys naics_2d: markupest mkup_hall, inputs(K L M S E) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) go(Y) deltago(deltaY) instruments(deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP)}
{phang2}{it:({stata "bys naics_2d: markupest mkup_hall, inputs( K L M S E ) deltavars( deltaK deltaL deltaM deltaS deltaE) method(hall) go(Y) deltago(deltaY) instruments( deltainstEq deltainstRD deltainstSh deltainstSo deltainstOP)":click to run})}

{pstd}Roeger - time invariant{p_end}
{tab}{cmd:. bys naics_2d: markupest mkup_roeg, inputs(L M S E  K) method(roeger) go(Y) }
{phang2}{it:({stata "bys naics_2d: markupest mkup_roeg, inputs(L M S E  K) method(roeger) go(Y)":click to run})}

{pstd}Collapse to compare{p_end}
{tab}{cmd:. collapse (mean) mkup_h mkup_roeg, by(naics_2d)}
{phang2}{it:({stata "collapse (mean) mkup_h mkup_roeg, by(naics_2d)":click to run})}

{title:References}

{phang}
Rovigatti, G., 2020, 
Markup estimation using Stata: the Micro and the Macro approaches with markupest.
Working Paper

{phang}
Rovigatti, G. and Mollisi, V., 2018, 
Theory and practice of total-factor productivity estimation: The control function approach using Stata.
The Stata Journal, 18, (3), pp.618-662

{phang}
De Loecker, J., and Warzynski, F., 2012, 
Markups and firm-level export status. 
American economic review 102(6), pp.2437–71

{phang}
Hall, R. E., 1988, 
The relation between price and marginal cost in US industry. 
Journal of political Economy 96(5), pp.921–947

{phang}
Hall, R. E., 2018, 
New evidence on the markup of prices over marginal costs and the role of mega-firms in the us economy. 
NBER Working Paper

{phang}
Roeger, W., 1995, 
Can imperfect competition explain the difference between primal and dual productivity measures? Estimates for US manufacturing 
Journal of political Economy 103(2), pp.316–330

{title:Saved results}

{pstd}
{cmd:markupest} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: General}{p_end}
{synopt:{cmd:e(cmd)}}markupest{p_end}
{synopt:{cmd:e(markupvar)}}name of estimated markups variable{p_end}
{synopt:{cmd:e(markuptype)}}type of estimated markups - micro or macro{p_end}
{synopt:{cmd:e(method)}}estimation method: dlw, hall, or roeger{p_end}
{synopt:{cmd:e(id)}}panel varname{p_end}
{synopt:{cmd:e(t)}}time varname{p_end}
{synopt:{cmd:e(w)}}regression weights{p_end}

{p2col 5 20 24 2: Micro}{p_end}
{synopt:{cmd:e(inputvar)}}varname of the variable input used for markup estimation(s){p_end}
{synopt:{cmd:e(output)}}varname of the output variable(s){p_end}
{synopt:{cmd:e(free)}}varname of the free variable(s){p_end}
{synopt:{cmd:e(state)}}varname of the state variable(s){p_end}
{synopt:{cmd:e(proxy)}}varname of the proxy variable(s){p_end}
{synopt:{cmd:e(PFest_method)}}PF estimation method (op, lp, wrdg, mr){p_end}
{synopt:{cmd:e(prodestOPTS)}}string with all prodest supplementary options{p_end}
{synopt:{cmd:e(corrected)}}markups corrected for first-stage residuals{p_end}

{p2col 5 20 24 2: Macro}{p_end}
{synopt:{cmd:e(go)}}varname of the gross output variable{p_end}
{synopt:{cmd:e(inputs)}}varnames of input variables{p_end}
{synopt:{cmd:e(deltago)}}varname of gross output changes. Hall method only{p_end}
{synopt:{cmd:e(deltavars)}}varnames of the changes in input variables. Hall method only{p_end}
{synopt:{cmd:e(pgo)}}varname of changes in gross output price. Hall method only{p_end}
{synopt:{cmd:e(pricevars)}}varnames of changes in input variable prices. Hall method only{p_end}
{synopt:{cmd:e(instruments)}}varname of the instruments. Hall method only{p_end}
{synopt:{cmd:e(timevarying)}}time-varying version of the Hall method{p_end}

{title:Authors}

{pstd}Gabriele Rovigatti{p_end}
{pstd}Bank of Italy{p_end}
{pstd}Rome, Italy{p_end}
{pstd}gabriele.rovigatti@gmail.com{p_end}

{title:Also see}

{psee}
Online: {helpb prodest}, {helpb prodest_p}{p_end}

{phang}
Rovigatti, G. and Mollisi, V., 2018, 
Theory and practice of total-factor productivity estimation: The control function approach using Stata.
The Stata Journal, 18, (3), pp.618-662

{title:Bug Reporting}

{psee}
{opt markupest} is part of an ongoing project. Please submit bugs, comments and suggestions via email to:

{pstd}	gabriele.rovigatti@gmail.com

{psee}
Please try to include steps to reproduce issues and possibly the output of the Stata command -creturn list-.


