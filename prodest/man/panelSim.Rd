\name{panelSim}
\alias{panelSim}
\title{
  Simulate Panel dataset
}
%%%%%%%%% DESCRIPTION %%%%%%%%%%%
\description{
  \code{panelSim()} produces a N*T balanced panel dataset of firms' production. In particular, it returns a \code{data.frame} with free, state and proxy variables aimed at performing Monte Carlo simulations on productivity-related models.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  panelSim(N = 1000, T = 100, alphaL = .6, alphaK = .4, DGP = 1,
           rho = .7, sigeps = .1, sigomg = .3, rholnw = .3, seed = 123456)
}
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
\arguments{
  \item{N}{
    the number of firms. By default \code{N=1000}}
  %
  \item{T}{
    the total time span to be simulated. Only a fraction (the last 10\% of observations) will be returned. By default \code{T=100}}
  %
  \item{alphaL}{
    the parameter of the free variable. By default \code{alphaL=.6}}
  %
  \item{alphaK}{
    the parameter of the state variable. By default \code{alphaK=.4}}
  %
  \item{DGP}{
    Type of DGP; accepts 1, 2 or 3. They differ in terms of shock to wages (0 or 0.1), \eqn{\Delta} (0 or 0.5) and shock to labor (0 or 0.37). See \eqn{details}. By default \code{DGP=1}.}
  %
  \item{rho}{
    the AR(1) coefficient for omega. By default \code{rho=0.7}}
  %
  \item{sigeps}{
    the standard deviation of epsilon. See \eqn{details}. By default \code{sigeps = .1}.}
  %
  \item{sigomg}{
    the standard deviation of the innovation to productivity \eqn{\omega}. By default \code{sigomg = .3}.}
  %
  \item{rholnw}{
    AR(1) coefficient for log(wage). By default \code{rholnw=.3}.}
  %
  \item{seed}{
    seed set when the routine starts. By default \code{seed = 123456}.}
  %
}
%%%%%%%%% DETAILS %%%%%%%%%%%
\details{
  \code{panelSim()} is the R implementation of the DGP written by Ackerberg, Caves and Frazer (2015).
}
%%%%%%%%% VALUE %%%%%%%%%%%
\value{
  \eqn{panelSim()} returns a \code{data.frame} with 7 variables:
  \itemize{
    \item{\eqn{idvar}}{
      ID codes from 1 to N (by default \eqn{N = 1000}).
    }
    \item{\eqn{timevar}}{
      time variable ranging 1 to \eqn{round(T*0.1)} (by default \eqn{T = 100} and \eqn{max(timevar) = 10}).
    }
    \item{\eqn{Y}}{
      log output value added variable
    }
    \item{\eqn{sX}}{
      log state variable
    }
    \item{\eqn{fX}}{
      log free variable
    }
    \item{\eqn{pX1}}{
      log proxy variable - no measurement error
    }
    \item{\eqn{pX2}}{
      log proxy variable - \eqn{\sigma_{measurementerror}= .1}
    }
    \item{\eqn{pX3}}{
      log proxy variable - \eqn{\sigma_{measurementerror}= .2}
    }
    \item{\eqn{pX4}}{
      log proxy variable - \eqn{\sigma_{measurementerror}= .5}
    }
  }
}
%%%%%%%%% AUTHOR %%%%%%%%%%%
\author{Gabriele Rovigatti}
%%%%%%%%% REFERENCES %%%%%%%%%%%
\references{
  Ackerberg, D., Caves, K. and Frazer, G. (2015).
  "Identication properties of recent production function estimators."
  Econometrica, 83(6), 2411-2451.
}
%%%%%%%%% EXAMPLES %%%%%%%%%%%
\examples{

  require(prodest)

  ## Simulate a dataset with 1000 firms (T = 100). \code{Panelsim()} delivers the last 10\% of usable time per panel.

  panel.data <- panelSim()
  attach(panel.data)

  ## Estimate various models
  ACF.fit <- prodestACF(Y, fX, sX, pX2, idvar, timevar, theta0 = c(.5,.5))
  LP.fit <- prodestLP(Y, fX, sX, pX2, idvar, timevar)
  WRDG.fit <- prodestWRDG(Y, fX, sX, pX3, idvar, timevar, R = 5)
  
  ## print results in lateX tabular format 
  printProd(list(LP.fit, ACF.fit, WRDG.fit))
}
