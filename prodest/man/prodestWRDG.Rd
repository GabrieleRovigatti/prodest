\name{prodestWRDG}
\alias{prodestWRDG}
\title{
  Estimate productivity - Wooldridge method
}
\description{
  The \code{prodestWRDG()} function accepts at least 6 objects (id, time, output, free, state and proxy variables), and returns a \code{prod} object of class \code{S4} with three elements: (i) a list of model-related objects, (ii) a list with the data used in the estimation and estimated vectors of first-stage residuals, and (iii) a list with the estimated parameters and their bootstrapped standard errors.
}
\usage{
  prodestWRDG(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL, seed = 123456,
                tol = 1e-100, theta0 = NULL, cluster = NULL)
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
    and \eqn{u_{it}} is a random shock component assumed to be uncorrelated with the technical efficiency, the state variables in \eqn{k_{it}} and the lagged free variables \eqn{w_{it-1}}.
    Wooldridge method allows to jointly estimate OP/LP two stages jointly in a system of two equations. It relies on the following set of assumptions:
      \itemize{
        \item a) \eqn{\omega_{it} = g(x_{it} , p_{it})}: productivity is an unknown function \eqn{g(.)} of state and a proxy variables;
\item b) \eqn{E(\omega_{it} | \omega_{it-1)}=f[\omega_{it-1}]}, productivity is an unknown function \eqn{f[.]} of lagged productivity, \eqn{\omega_{it-1}}.
      }
   Under the above set of assumptions, It is possible to construct a system gmm using the  vector of residuals from
   \itemize{
    \item \eqn{r_{1it} = y_{it} - alpha - w_{it}\beta - x_{it}\gamma - g(x_{it} , p_{it}) }
    \item \eqn{r_{2it} = y_{it} - alpha - w_{it}\beta - x_{it}\gamma - f[g(x_{it-1} , p_{it-1})]}
   }
  where the unknown function \eqn{f(.)} is approximeted by a n-th order polynomial and \eqn{g(x_{it} , m_{it}) = \lambda_0 + c(x_{it} , m_{it})\lambda}. In particular, \eqn{g(x_{it} , m_{it})} is a linear combination of functions in \eqn{(x_{it} , m_{it})}
and \eqn{c_{it}} are the addends of this linear combination. The residuals eqn{r_{it}} are used to set the moment conditions
  \itemize{
    \item \eqn{E(Z_{it}*r_{it}) =0}
  }
  with the following set of instruments:
  \itemize{
    \item \eqn{Z1_{it} = (1, w_{it}, x_{it}, c_{it})}
    \item \eqn{Z2_{it} = (w_{it-1}, c_{it}, c_{it})}
  }
  }
%%%%%%%%% VALUE %%%%%%%%%%%
  \value{
The output of the function \code{prodestWRDG} is a member of the \code{S3} class \pkg{prod}. More precisely, is a list (of length 3) containing the following elements:

    \code{Model}, a list containing:
    \itemize{
      \item \code{method:} a string describing the method ('WRDG').
      \item \code{elapsed.time:} time elapsed during the estimation.
      \item \code{seed:} the seed set at the beginning of the estimation.
      \item \code{opt.outcome:} optimization outcome.
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
    }

    \code{Estimates}, a list containing:
    \itemize{
      \item \code{pars:} the vector of estimated coefficients.
      \item \code{std.errors:} the vector of bootstrapped standard errors.
    }

    Members of class \code{prod} have an \code{omega} method returning a numeric object with the estimated productivity - that is: \eqn{\omega_{it} = y_{it} - (\alpha + w_{it}\beta + k_{it}\gamma)}.
    \code{FSres} method returns a numeric object with the residuals of the first stage regression, while \code{summary}, \code{show} and \code{coef} methods are implemented and work as usual.
  }
%%%%%%%%% NOTE %%%%%%%%%%%
%  \note{
%    \code{prodestWRDG()} proceeds with the estimation
%  }
%%%%%%%%% AUTHOR %%%%%%%%%%%
  \author{
    Gabriele Rovigatti
  }
%%%%%%%%% REFERENCES %%%%%%%%%%%
  \references{
    Wooldridge, J M (2009).
    "On estimating firm-level production functions using proxy variables to control for unobservables."
    Economics Letters, 104, 112-114.
  }
%%%%%%%%% EXAMPLES %%%%%%%%%%%
  \examples{
    data("chilean")

    # we fit a model with two free (skilled and unskilled), one state (capital) and one proxy variable (electricity)

    WRDG.fit <- prodestWRDG(d$Y, fX = cbind(d$fX1, d$fX2), d$sX, d$pX, d$idvar, d$timevar)

    # show results
    WRDG.fit

    \dontrun{
      # estimate a panel dataset - DGP1, various measurement errors - and run the estimation
      sim <- panelSim()

      WRDG.sim1 <- prodestWRDG(sim$Y, sim$fX, sim$sX, sim$pX1, sim$idvar, sim$timevar)
      WRDG.sim2 <- prodestWRDG(sim$Y, sim$fX, sim$sX, sim$pX2, sim$idvar, sim$timevar)
      WRDG.sim3 <- prodestWRDG(sim$Y, sim$fX, sim$sX, sim$pX3, sim$idvar, sim$timevar)
      WRDG.sim4 <- prodestWRDG(sim$Y, sim$fX, sim$sX, sim$pX4, sim$idvar, sim$timevar)

      # show results in .tex tabular format
      printProd(list(WRDG.sim1, WRDG.sim2, WRDG.sim3, WRDG.sim4), parnames = c('Free','State'))
    }
  }

