\name{chilean}
\docType{data}
\alias{chilean}
\title{Data: Chilean firm-level production data 1986-1996}
\description{
  Sectoral subsample of Chilean firm-level production data 1986-1996.
}
\usage{data("chilean")}
\format{A \link{data.frame} object containing 9 variables with production-related data.}
%%%%%%%%% VALUEs %%%%%%%%%%%
\value{
  \item{Y}{
    vector of log(outcome) - Value added.}
  %
  \item{sX}{
    vector of log(capital).}
  %
  \item{fX}{
    matrix of log(skilled labor) and log(unskilled labor).}
  %
  \item{cX}{
    vector of log(water).}
  %
  \item{pX}{
    vector of log(electricity).}
  %
  \item{inv}{
    vector of log(investment).}
  %
  \item{idvar}{
    vector of panel identifier.}
  %
  \item{timevar}{
    vector of time.}
  %
}
%%%%%%%%% REFERENCES %%%%%%%%%
\references{
  \url{http://www.ine.cl/canales/chile_estadistico/estadisticas_economicas/industria/series_estadisticas/series_estadisticas_enia.php}
}
%%%%%%%%% KEYWORDS %%%%%%%%%
\keyword{datasets}
