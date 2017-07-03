\name{prodestOP}
\alias{prodestOP}
\title{
	Estimate productivity - Olley-Pakes method
}
\description{
The \code{prodestOP()} function accepts at least 6 objects (id, time, output, free, state and proxy variables), and returns a \code{prod} object of class \code{S4} with three elements: (i) a list of model-related objects, (ii) a list with the data used in the estimation and estimated vectors of first-stage residuals, and (iii) a list with the estimated parameters and their bootstrapped standard errors .
}
\usage{
  prodestOP(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL,
            opt = 'optim', theta0 = NULL, seed = 123456, cluster = NULL, tol = 1e-100)
}
\arguments{
  \item{Y      }{
  the vector of value added log output.}
  %
  \item{fX     }{
  the vector/matrix/dataframe of log free variables.}
  %
  \item{sX     }{
  the vector/matrix/dataframe of log state variables.}
  %
  \item{pX    }{
  the vector/matrix/dataframe of log proxy variables.}
  %
  \item{cX     }{
  the vector/matrix/dataframe of control variables. By default \code{cX= NULL}.}
  %
  \item{idvar  }{
  the vector/matrix/dataframe identifying individual panels.}
  %
  \item{timevar}{
  the vector/matrix/dataframe identifying time.}
  %
  \item{R      }{
  the number of block bootstrap repetitions to be performed in the standard error estimation. By default \code{R = 20}.}
  %
  \item{opt    }{
  a string with the optimization algorithm to be used during the estimation. By default \code{opt = 'optim'}.}
  %
  \item{theta0 }{
  a vector with the second stage optimization starting points. By default \code{theta0 = NULL} and the optimization is run starting from the first stage estimated parameters + \eqn{N(\mu=0,\sigma=0.01)} noise.}
  %
  \item{cluster}{
    an object of class \code{"SOCKcluster"} or \code{"cluster"}. By default \code{cluster = NULL}.}
  %
  \item{seed}{
    seed set when the routine starts. By default \code{seed = 123456}.}
  %
  \item{tol}{
    optimizer tolerance. By default \code{tol = 1e-100}.}
  %
}
%%%%%%%%% DETAILS %%%%%%%%%%%
\details{
  Consider a Cobb-Douglas production technology for firm \eqn{i} at time \eqn{t}
  \itemize{
    \item \eqn{y_{it} = \alpha + w_{it}\beta + k_{it}\gamma + \omega_{it} + \epsilon_{it}}
  }
  where \eqn{y_{it}} is the (log) output, w_{it} a 1xJ vector of (log) free variables, k_{it} is a 1xK vector of state variables and \eqn{\epsilon_{it}} is a normally distributed idiosyncratic error term.
  The unobserved technical efficiency parameter \eqn{\omega_{it}} evolves according to a first-order Markov process:
  \itemize{
    \item \eqn{\omega_{it} = E(\omega_{it} | \omega_{it-1}) + u_{it} = g(\omega_{it-1}) + u_{it}}
  }
  and \eqn{u_{it}} is a random shock component assumed to be uncorrelated with the technicalefficiency, the state variables in \eqn{k_{it}} and the lagged free variables \eqn{w_{it-1}}.
  The OP method relies on the following set of assumptions:
  \itemize{
    \item a) \eqn{i_{it} = i(k_{it},\omega_{it})} - investments are a function of both the state variable and the technical efficiency parameter;
    \item b) \eqn{i_{it}} is strictly monotone in \eqn{\omega_{it}};
    \item c) \eqn{\omega_{it}} is scalar unobservable in \eqn{i_{it} = i(.)} ;
    \item d) the levels of \eqn{i_{it}} and \eqn{k_{it}} are decided at time \eqn{t-1}; the level of the free variable, \eqn{w_{it}}, is decided after the shock \eqn{u_{it}} realizes.
  }
  Assumptions a)-d) ensure the invertibility of \eqn{i_{it}} in \eqn{\omega_{it}} and lead to the partially identified model:
  \itemize{
    \item \eqn{y_{it} = \alpha + w_{it}\beta + k_{it}\gamma + h(i_{it}, k_{it}) + \epsilon_{it} = \alpha + w_{it}\beta  + \phi(i_{it}, k_{it}) + \epsilon_{it} }
  }
  which is estimated by a non-parametric approach - First Stage.
  Exploiting the Markovian nature of the productivity process one can use assumption d) in order to set up the relevant moment conditions and estimate the production function parameters - Second stage.
  Exploting the resisual \eqn{e_{it}} of:
  \itemize{
    \item \eqn{y_{it} - w_{it}\hat{\beta} = \alpha + k_{it}\gamma + g(\omega_{it-1}, \chi_{it}) + \epsilon_{it} }
  }
  and \eqn{g(.)} is typically left unspecified and approximated by a \eqn{n^{th}} order polynomial and \eqn{\chi_{it}} is an inidicator function for the attrition in the market.
}
%%%%%%%%% VALUE %%%%%%%%%%%
\value{
    The output of the function \code{prodestOP} is a member of the \code{S3} class \pkg{prod}. More precisely, is a list (of length 3) containing the following elements:

    \code{Model}, a list containing:
    \itemize{
      \item \code{method:} a string describing the method ('OP').
      \item \code{boot.repetitions:} the number of bootstrap repetitions used for standard errors' computation.
      \item \code{elapsed.time:} time elapsed during the estimation.
      \item \code{theta0:} numeric object with the optimization starting points - second stage.
      \item \code{opt:} string with the optimization routine used - 'optim', 'solnp' or 'DEoptim'.
      \item \code{seed:} the seed set at the beginning of the estimation.
      \item \code{opt.outcome:} optimization outcome.
      \item \code{FSbetas:} first stage estimated parameters.
    }

    \code{Data}, a list containing:
    \itemize{
      \item \code{Y:} the vector of value added log output.
      \item \code{free:} the vector/matrix/dataframe of log free variables.
      \item \code{state:} the vector/matrix/dataframe of log state variables.
      \item \code{proxy:} the vector/matrix/dataframe of log proxy variables.
      \item \code{control:} the vector/matrix/dataframe of log control variables.
      \item \code{idvar:} the vector/matrix/dataframe identifying individual panels.
      \item \code{timevar:} the vector/matrix/dataframe identifying time.
      \item \code{FSresiduals:} numeric object with the residuals of the first stage.
    }

    \code{Estimates}, a list containing:
    \itemize{
      \item \code{pars:} the vector of estimated coefficients.
      \item \code{std.errors:} the vector of bootstrapped standard errors.
    }

    Members of class \code{prod} have an \code{omega} method returning a numeric object with the estimated productivity - that is: \eqn{\omega_{it} = y_{it} - (\alpha + w_{it}\beta + k_{it}\gamma)}.
    \code{FSres} method returns a numeric object with the residuals of the first stage regression, while \code{summary}, \code{show} and \code{coef} methods are implemented and work as usual.
}
%%%%%%%%% AUTHOR %%%%%%%%%%%
\author{Gabriele Rovigatti}
%%%%%%%%% REFERENCES %%%%%%%%%%%
\references{
  Olley, S G and Pakes, A (1996).
  "The dynamics of productivity in the telecommunications equipment industry."
  Econometrica, 64(6), 1263-1297.
}
%%%%%%%%% EXAMPLES %%%%%%%%%%%
\examples{
    require(prodest)

    ## Chilean data on production.The full version is Publicly available at
    ## http://www.ine.cl/canales/chile_estadistico/estadisticas_economicas/industria/series_estadisticas/series_estadisticas_enia.php

    data(chilean)

    # we fit a model with two free (skilled and unskilled), one state (capital) and one proxy variable (electricity)

    OP.fit <- prodestOP(d$Y, fX = cbind(d$fX1, d$fX2), d$sX, d$inv, d$idvar, d$timevar)
    OP.fit.solnp <- prodestOP(d$Y, fX = cbind(d$fX1, d$fX2), d$sX, d$inv, d$idvar, d$timevar, opt='solnp')
    OP.fit.control <- prodestOP(d$Y, fX = cbind(d$fX1, d$fX2), d$sX, d$inv, d$idvar, d$timevar, cX = d$cX)

    # show results
    summary(OP.fit)
    summary(OP.fit.solnp)
    summary(OP.fit.control)

    # show results in .tex tabular format
     printProd(list(OP.fit, OP.fit.solnp, OP.fit.control))
}
