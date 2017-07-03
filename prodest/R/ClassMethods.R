# set the "prod" object class #
setClass("prod", representation(Model = "list", Data = "list", Estimates = "list"))
# end of object classs setting #

# coefficients #
setMethod("coef", signature(object = "prod"), function(object, iBurnPeriod = NULL ) object@Estimates$pars)
# end coefficients #

# Report the Omega estimate #
setGeneric('omega', function(object, ...) {object}) # set the generic "omega" function for "prod" objects
setMethod("omega", signature(object = "prod"), function(object, iBurnPeriod = NULL ){
  Omega <- object@Data$Y - cbind(object@Data$free, object@Data$state, object@Data$control) %*% object@Estimates$pars
  return(Omega)
})
# end of omega #

# First stage residuals #
setGeneric('FSres', function(object, ...) {object}) # set the generic "FSres" function for "prod" objects
setMethod("FSres", signature(object = "prod"), function(object, iBurnPeriod = NULL ) object@Data$FSresiduals)
# end FS residuals #

# Show method #
setMethod('show', signature(object = 'prod'), function(object){

  Pars = object@Estimates$pars
  nPars = names(Pars)
  st.err = object@Estimates$std.errors
  Method = object@Model$method

  cat("\n-------------------------------------------------------")
  cat("\n-            Production Function Estimation           -")
  cat("\n-------------------------------------------------------")
  cat(paste("\n                   Method:   ", Method, "             "))
  cat("\n-------------------------------------------------------")
  cat(paste("\n                      ", paste("  ", nPars, collapse = "    ", " ", sep = "")))
  cat(paste("\nEstimated Parameters: ", paste(" ", round(Pars, digits = 3), collapse = "   ", " ", sep = "")))
  cat(paste("\n                      ", paste("(", round(st.err, digits = 3), collapse = "   ", ")", sep = "")))
  cat("\n-------------------------------------------------------")
})
# end of show method

# Show method #
setMethod('summary', signature(object = 'prod'), function(object){

  Pars = object@Estimates$pars
  namePars = names(Pars)
  st.err = object@Estimates$std.errors
  Method = object@Model$method
  N = nrow(object@Data$Y)
  Time = object@Model$elapsed.time
  nCores = object@Model$nCores
  opt = object@Model$opt
  R = object@Model$boot.repetitions
  theta0 = object@Model$theta0

  cat("\n-------------------------------------------------------------")
  cat("\n-               Production Function Estimation              -")
  cat("\n-------------------------------------------------------------")
  cat(paste("\n                   Method :   ", Method, "             "))
  cat("\n-------------------------------------------------------------")
  cat(paste("\n                            ", paste("  ", namePars, collapse = "    ", " ", sep = "")))
  cat(paste("\nEstimated Parameters      : ", paste(" ", round(Pars, digits = 3), collapse = "   ", " ", sep = "")))
  cat(paste("\n                            ", paste("(", round(st.err, digits = 3), collapse = "   ", ")", sep = "")))
  cat("\n-------------------------------------------------------------")
  cat(paste("\nN                         :  ", N, sep = "") )
  cat("\n-------------------------------------------------------------")
  if (!(Method == 'WRDG')){
    cat(paste("\nBootstrap repetitions     : ", R))
    FS.betas = object@Model$FSbetas[2 : (1 + length(Pars))]
    cat(paste("\n1st Stage Parameters      :", paste(" ", round(FS.betas, digits = 3), collapse = "   ", " ", sep = "")))
    cat(paste("\nOptimizer                 : ", opt))
    if (!is.null(theta0)){
      cat(paste("\n2nd Stage Start Points    :", paste(" ", theta0, collapse = "   ", " ", sep = "")))
    }
  }
  cat("\n-------------------------------------------------------------")
  cat(paste("\nElapsed Time              :  ", round(as.double(Time, units = 'mins'), digits = 2), " mins", sep = "") )
  if (!is.null(nCores)){
    cat(paste("\n# Cores                 : ", nCores, sep = "") )
  }
  cat("\n-------------------------------------------------------------")
})
# end of show method
