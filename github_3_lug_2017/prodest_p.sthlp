{smcl}
{* *! version 1.0.1  06Jun2017}{...}

{cmd:help prodest predict} {right:also see: {helpb prodest}  }
{hline}

{title:Title}

{p 4 16 2}
{cmd:prodest predict} {hline 2} Postestimation tool for prodest{p_end}


{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} [{newvar}] {ifin} [{cmd:,} {cmdab:resid:uals} {cmdab:exp:onential} {cmdab:par:ameters}]

{synoptset 28 tabbed}{...}
{synopthdr :options}
{synoptline}
{synopt :{opt resid:uals}}stores the residuals - log(TFP) - in {newvar}{p_end}
{synopt :{opt exp:onential}}stores the exponentiated residuals - TFP - in {newvar}{p_end}
{synopt :{opt par:ameters}}returns the estimated input elasticities{p_end}
{synoptline}
{p2colreset}{...}

{title:Options for predict}

{dlgtab:Cobb-Douglas PF}

{phang}
{opt residuals:} log(TFP) values calculated from the log production function as:
omega_it = y_it - beta_{l} * l - beta_{k} * k

{phang}
{opt exponentiated:} TFP values calculated as the exponential of the residuals from the log production function as:
exp(omega_it) = exp(y_it - beta_{l} * l - beta_{k} * k)

{phang}
{opt parameters:} The estimated parameters for free, state and control variables.

{dlgtab:Translog PF}

{phang}
{opt residuals:} log(TFP) values calculated from the log production function as:
omega_it = y_it - beta_{l} * l - beta_{k} * k - beta_{ll} * l^2 - beta_{kk} * k^2 - beta_{lk} * (l * k)

{phang}
{opt exponentiated:} TFP values calculated as the exponential of the residuals from the log production function as:
exp(omega_it) = exp(y_it - beta_{l} * l - beta_{k} * k - beta_{ll} * l^2 - beta_{kk} * k^2 - beta_{lk} * (l * k))

{phang}
{opt parameters:} The estimated elasticities for free and state variables. With translog PF, the elasticities (beta) are defined as:
hat(beta)_{l} = E( beta_{l} + 2 beta_{ll} * l + beta_{lk} * k )
hat(beta)_{k} = E( beta_{k} + 2 beta_{kk} * k + beta_{lk} * l )

{marker examples}{...}
{title:Examples}
{tab}{cmd:. insheet using "https://raw.githubusercontent.com/GabBrock/prodest/master/prodest.csv", names clear}}
{phang2}{it:({stata "insheet using https://raw.githubusercontent.com/GabBrock/prodest/master/prodest.csv, names clear":load data})}

{tab}{cmd:. xtset id year, y /* not run */}

{pstd}Levinsohn and Petrin method{p_end}
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) opt(dfp) reps(50) id(id) t(year) fsresiduals(fs_lp)}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) opt(dfp) reps(50) id(id) t(year) fsres(sf_lp)":click to run})}

{tab}{cmd:. predict lpfit, residuals}
{phang2}{it:({stata "predict lpfit, residuals":click to run})}

{tab}Ackerberg, Caves and Frazer with translog
{tab}{cmd:. prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) acf reps(20) id(id) t(year) translog}
{phang2}{it:({stata "prodest log_y, free(log_lab1 log_lab2) state(log_k) proxy(log_materials) va met(lp) acf reps(20) id(id) t(year) translog":click to run})}

{tab}{cmd:. predict, parameters}
{phang2}{it:({stata "predict, parameters":click to run})}

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

