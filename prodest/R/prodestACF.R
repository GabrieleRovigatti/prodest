############ ACKERBERG-CAVES-FRAZER ###############

# function to estimate ACF model #
prodestACF <- function(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL, opt = 'optim',
                       theta0 = NULL, seed = 123456, cluster = NULL){
  set.seed(seed)
  Start = Sys.time() # start tracking time
  Y <- checkM(Y) # change all input to matrix
  fX <- checkM(fX)
  sX <- checkM(sX)
  pX <- checkM(pX)
  idvar <- checkM(idvar)
  timevar <- checkM(timevar)
  snum <- ncol(sX) # find the number of input variables
  fnum <- ncol(fX)
  if (!is.null(cX)) {cX <- checkM(cX); cnum <- ncol(cX)} else {cnum <- 0} # if is there any control, take it into account, else fix the number of controls to 0
  if (length(theta0) != cnum + fnum + snum & !is.null(theta0)){
    stop(paste0('theta0 length (', length(theta0), ') is inconsistent with the number of parameters (', cnum + fnum + snum, ')'), sep = '')
  }
  polyframe <- data.frame(fX,sX,pX) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(mod, fX^2, sX^2, pX^2) # generate a polynomial of the desired level
  lag.sX = sX # generate sX lags
  for (i in 1:snum) {
    lag.sX[, i] = lagPanel(sX[, i], idvar = idvar, timevar = timevar)
  }
  lag.fX = fX # generate fX lags
  for (i in 1:fnum) {
    lag.fX[, i] = lagPanel(fX[, i], idvar = idvar, timevar = timevar)
  }
  if (!is.null(cX)) { # generate the matrix of data
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX, sX),
                                 Xt = data.frame(fX, sX), lX = data.frame(lag.fX, lag.sX),
                                 cX = data.frame(cX), regvars = regvars))
  } else {
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX, sX),
                                 Xt = data.frame(fX, sX), lX = data.frame(lag.fX,lag.sX), regvars = regvars))
  }
  betas <- finalACF(ind = TRUE, data = data, fnum = fnum, snum = snum, cnum = cnum, opt = opt, theta0 = theta0)
  # if (betas$opt.outcome$convergence != 0){
  #   warning('Second Stage convergence not achieved')
  # }
  boot.indices <- block.boot.resample(idvar, R) # generates a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    nCores = NULL
    boot.betas <- matrix(unlist(
      lapply(boot.indices, finalACF, data = data, fnum = fnum, snum = snum, cnum = cnum, opt = opt,
             theta0 = theta0, boot = TRUE)), ncol = fnum + snum + cnum, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  } else {
    nCores = length(cluster)
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix( unlist( parLapply(cl = cluster, boot.indices, finalACF, data = data, fnum = fnum, snum = snum,
                                            cnum = cnum, opt = opt, theta0 = theta0, boot = TRUE) ),
                          ncol = fnum + snum + cnum, byrow = TRUE ) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas, 2, sd, na.rm = TRUE) # calculate standard deviations
  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(cX)) {
    res.names <- c(res.names, colnames(cX, do.NULL = FALSE, prefix = 'cX'))
  }
  names(betas$betas) <- res.names # change results' names
  names(boot.errors) <- res.names # change results' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'ACF', FSbetas = NA, boot.repetitions = R, elapsed.time = elapsedTime, theta0 = theta0,
                          opt = opt, seed = seed, opt.outcome = betas$opt.outcome, nCores = nCores),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar,
                         FSresiduals = betas$FSresiduals),
             Estimates = list(pars = betas$betas, std.errors = boot.errors))
  return(out)
}
# end of prodestACF #

# function to estimate and bootstrap ACF #
finalACF <- function(ind, data, fnum, snum, cnum, opt, theta0, boot = FALSE){
  if (sum(as.numeric(ind)) == length(ind)){ # if the ind variable is not always TRUE
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  data <- data[ind,] # change the index according to bootstrapped indices
  first.stage <- lm(data[,'Y', drop = FALSE] ~ data[, grepl('regvars', colnames(data)), drop = FALSE], na.action = na.exclude)
  phi <- fitted(first.stage) # generate the fitted values of the first stage
  if (is.null(theta0)) {
    theta0 <- coef(first.stage)[2:(1 + snum + fnum + cnum)] + rnorm((snum + fnum), 0, 0.01)
  } # use the first stage + noise results as starting points in case the user did not specify other
  newtime <- data[,'timevar', drop = FALSE]
  rownames(phi) <- NULL
  rownames(newtime) <- NULL
  lag.phi <- lagPanel(idvar = newid, timevar = newtime, value = phi) #  # lag fitted values
  Z <- data[, grepl('Z', colnames(data)), drop = FALSE]
  X <- data[, grepl('Xt', colnames(data)), drop = FALSE]
  lX <- data[, grepl('lX', colnames(data)), drop = FALSE]
  tmp.data <- model.frame(Z ~ X + lX + phi + lag.phi)
  W <- solve(crossprod(tmp.data$Z)) / nrow(tmp.data$Z)
  if (opt == 'optim'){
    try.out <- try(optim(theta0, gACF, method = "BFGS", mZ = tmp.data$Z, mW = W, mX = tmp.data$X,
                         mlX = tmp.data$lX,
                         vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$par
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA,(snum + fnum), 1)
      opt.outcome <- list(convergence = 999)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'DEoptim'){
    try.out <- try(DEoptim(gACF, lower = theta0, upper = rep.int(1,length(theta0)), mZ = tmp.data$Z,
                           mW = W, mX = tmp.data$X, mlX = tmp.data$lX,
                           vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi,
                           control = DEoptim.control(trace = FALSE)), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$optim$bestmem
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA, (snum + fnum), 1)
      opt.outcome <- list(convergence = 99)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'solnp'){
    try.out <- try(solnp(theta0, gACF, mZ = tmp.data$Z, mW = W, mX = tmp.data$X, mlX = tmp.data$lX,
                         vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi,
                         control = list(trace = FALSE)), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$pars
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA,(snum + fnum), 1)
      opt.outcome <- list(convergence = 999)
    } # error handling: if the optimization fails return missing values
  }
  if (boot == FALSE){ # for the baseline estimation we want info on optimization, too
    return(list(betas = betas, opt.outcome = opt.outcome, FSresiduals = resid(first.stage)))
  } else {
    return(betas)
  }
}
# end of ACF final function #

# function to run the GMM estimation for ACF #
gACF <- function(theta, mZ, mW, mX, mlX, vphi, vlag.phi){
  Omega <- vphi - mX %*% theta
  Omega_lag <- vlag.phi - mlX %*% theta
  Omega_lag_pol <- cbind(1, Omega_lag, Omega_lag^2, Omega_lag^3)
  g_b <- solve(crossprod(Omega_lag_pol)) %*% t(Omega_lag_pol) %*% Omega
  XI <- Omega - Omega_lag_pol %*% g_b
  crit <- t(crossprod(mZ, XI)) %*% mW %*% (crossprod(mZ, XI))
  return(crit)
}
# end of GMM ACF #
