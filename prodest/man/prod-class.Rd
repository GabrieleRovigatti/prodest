\name{prod}
\Rdversion{0.1}
\docType{class}
\alias{prod}
\alias{prod-class}
\title{Class for Prodest Fitted object}
\description{
  Class for prodest fitted objects.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    %
    \item{\code{Model}:}{Object of class \code{list}. Contains information about the model and the optimization procedue:
        \itemize{
          \item \code{method}: \code{string} The method used in estimation.
          \item \code{FSbetas}: \code{numeric} First-stage estimated parameters.
          \item \code{boot.repetitions}: \code{numeric} Number of bootstrap repetitions.
          \item \code{elapsed.time}: \code{numeric} Time - in seconds - required for estimation.
          \item \code{theta0}: \code{numeric} Vector of Second-stage optimization starting points.
          \item \code{opt}: \code{string} Optimizer used for the Second-stage.
          \item \code{seed}: \code{numeric} seed set.
          \item \code{opt.outcome}: \code{list} Optimization outcome (depends on optimizer choice).
        }
    }
    %
    \item{\code{Data}:}{Object of class \code{list}. Contains:
      \itemize{
        \item \code{Y}: \code{numeric} Dependent variable - Value added.
        \item \code{free}: \code{matrix} Free variable(s).
        \item \code{state}: \code{matrix} State variable(s).
        \item \code{proxy}: \code{matrix} Proxy variable(s).
        \item \code{control}: \code{matrix} Control variable(s).
        \item \code{idvar}: \code{numeric} Panel identifiers.
        \item \code{timevar}: \code{numeric} Time identifiers.
        \item \code{FSresiduals}: \code{numeric} First-Stage residuals.
      }
    }
    %
    \item{\code{Estimates}:}{Object of class \code{list}. Contains:
      \itemize{
        \item \code{pars}: \code{numeric} Estimated parameters for the variables of interest.
        \item \code{std.errors}: \code{numeric} Estimated standard errors for the variables of interest.
      }
    }
  }
}
\section{Methods}{
  \itemize{
    \item \code{show} \code{signature(object = 'prod')}: Show table with the method, the estimated parameters and their standard errors.
    \item \code{summary} \code{signature(object = 'prod')}: Show table with method, parameters, std.errors and auxiliary information on model and optimization.
    \item \code{FSres} \code{signature(object = 'prod')}: Extract First-Stage residual vector.
    \item \code{omega} \code{signature(object = 'prod')}: Extract estimated productivity vector.
    \item \code{coef} \code{signature(object = 'prod')}: Extract estimated coefficients.
  }
}
\author{Gabriele Rovigatti}
\keyword{classes}
