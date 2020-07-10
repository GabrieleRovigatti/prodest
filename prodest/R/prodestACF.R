############ ACKERBERG-CAVES-FRAZER (WITH DE LOECKER, GOLDBERG, PAVCNIK AND KHANDELWAL INPUT PRICE CONTROLS) ###############

# function to estimate ACF model #
prodestACF <- function(Y, fX, sX, pX, idvar, timevar, zX = NULL, control = 'none', dum = F, G = 3, A = 3, R = 20, orth = F, opt = 'optim',
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
  if (!is.null(zX) & control == "2s") {zX <- checkM(zX); znum <- ncol(zX)} else {znum <- 0} # to determine number of variables later
  if (length(theta0) != znum + fnum + snum & !is.null(theta0)){
    stop(paste0('theta0 length (', length(theta0), ') is inconsistent with the number of parameters (', znum + fnum + snum, ')'), sep = '')
  }
  if (!is.null(zX)) {
    polyframe <- poly(matrix(cbind(fX,sX,pX), ncol = fnum + snum + ncol(pX)),degree=G,raw=!orth) # create (orthogonal / raw) polynomial of degree G
    regvars <- cbind(fX,sX,zX,pX,polyframe) # to make sure 1st degree variables come first (lm will drop the other ones from 1st-stage reg)
  } else { 
    polyframe <- poly(matrix(cbind(fX,sX,pX), ncol = fnum + snum + ncol(pX)),degree=G,raw=!orth)
    regvars <- cbind(fX,sX,pX,polyframe) # to make sure 1st degree variables come first (lm will drop the other ones from 1st-stage reg)
  } # create (orthogonal / raw) polynomial of degree G
  if (dum) {regvars <- cbind(regvars, factor(timevar))} # add time dummies to first stage
  lag.sX = sX # generate sX lags
  for (i in 1:snum) {
    lag.sX[, i] = lagPanel(sX[, i], idvar = idvar, timevar = timevar)
  }
  lag.fX = fX # generate fX lags
  for (i in 1:fnum) {
    lag.fX[, i] = lagPanel(fX[, i], idvar = idvar, timevar = timevar)
  }
  if (control == "2s") { # generate zX lags only if we include controls in production function 
  lag.zX = zX # generate fX lags
  for (i in 1:znum) {
    lag.zX[, i] = lagPanel(zX[, i], idvar = idvar, timevar = timevar)
  }
  }
  if (!is.null(zX) & control == "fs") { # generate the matrix of data for case where controls only appear in first stage, as in De Loecker and Warzynski (2012)
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX),
                                 Xt = data.frame(fX,sX), lX = data.frame(lag.fX,lag.sX),
                                 zX = data.frame(zX), regvars = regvars))
  } else if (!is.null(zX) & control == "2s") {
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX,lag.zX),
                                 Xt = data.frame(fX,sX,zX), lX = data.frame(lag.fX,lag.sX,lag.zX),
                                 zX = data.frame(zX), regvars = regvars))
  } else {
    data <- as.matrix(data.frame(Y = Y, idvar = idvar, timevar = timevar, Z = data.frame(lag.fX,sX),
                                 Xt = data.frame(fX,sX), lX = data.frame(lag.fX,lag.sX), regvars = regvars))
  }
  betas <- finalACF(ind = TRUE, data = data, fnum = fnum, snum = snum, znum = znum, opt = opt, theta0 = theta0, A=A)
  if (betas$opt.outcome$convergence != 0){
    warning('Second Stage convergence not achieved')
  }
  boot.indices <- block.boot.resample(idvar, R) # generates a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    nCores = NULL
    boot.betas <- matrix(unlist(
      lapply(boot.indices, finalACF, data = data, fnum = fnum, snum = snum, znum = znum, opt = opt,
             theta0 = theta0, boot = TRUE)), ncol = fnum + snum + znum, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  } else {
    nCores = length(cluster)
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix( unlist( parLapply(cl = cluster, boot.indices, finalACF, data = data, fnum = fnum, snum = snum,
                                            znum = znum, opt = opt, theta0 = theta0, boot = TRUE) ),
                          ncol = fnum + snum + znum, byrow = TRUE ) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas, 2, sd, na.rm = TRUE) # calculate standard deviations
  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(zX) & control == "2s") {
    res.names <- c(res.names, colnames(zX, do.NULL = FALSE, prefix = 'zX'))
  }
  names(betas$betas) <- res.names # change results' names
  names(boot.errors) <- res.names # change results' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'ACF', FSbetas = NA, boot.repetitions = R, elapsed.time = elapsedTime, theta0 = theta0,
                          opt = opt, seed = seed, opt.outcome = betas$opt.outcome, nCores = nCores),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = zX, idvar = idvar, timevar = timevar,
                         FSresiduals = betas$FSresiduals),
             Estimates = list(pars = betas$betas, std.errors = boot.errors))
  return(out)
}
# end of prodestACF #

# function to estimate and bootstrap ACF #
finalACF <- function(ind, data, fnum, snum, znum, opt, theta0, A, boot = FALSE){
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
    theta0 <- coef(first.stage)[2:(1 + snum + fnum + znum)] + rnorm((snum + fnum + znum), 0, 0.01)
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
                         vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi, A = A,
                         control = list(maxit = 1500) # to 'guarantee' 2nd stage convergence
                         ), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$par
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA,(snum + fnum + znum), 1)
      opt.outcome <- list(convergence = 999)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'DEoptim'){
    try.out <- try(DEoptim(gACF, lower = theta0, upper = rep.int(1,length(theta0)), mZ = tmp.data$Z,
                           mW = W, mX = tmp.data$X, mlX = tmp.data$lX,
                           vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi, A = A,
                           control = DEoptim.control(trace = FALSE)), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$optim$bestmem
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA, (snum + fnum + znum), 1)
      opt.outcome <- list(convergence = 99)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'solnp'){
    try.out <- try(solnp(theta0, gACF, mZ = tmp.data$Z, mW = W, mX = tmp.data$X, mlX = tmp.data$lX,
                         vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi, A = A,
                         control = list(trace = FALSE)), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      betas <- try.out$pars
      opt.outcome <- try.out
    } else {
      betas <- matrix(NA,(snum + fnum + znum), 1)
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
gACF <- function(theta, mZ, mW, mX, mlX, vphi, vlag.phi, A){
  Omega <- vphi - mX %*% theta
  Omega_lag <- vlag.phi - mlX %*% theta
  Omega_lag_pol <- poly(Omega_lag,degree=A,raw=T) # create polynomial in omega for given degree
  Omega_lag_pol <- cbind(1, Omega_lag_pol)
  g_b <- solve(crossprod(Omega_lag_pol)) %*% t(Omega_lag_pol) %*% Omega
  XI <- Omega - Omega_lag_pol %*% g_b
  crit <- t(crossprod(mZ, XI)) %*% mW %*% (crossprod(mZ, XI))
  return(crit)
}
# end of GMM ACF #
