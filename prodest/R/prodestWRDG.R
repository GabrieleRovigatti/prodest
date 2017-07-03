############ WOOLDRIDGE ###############

# function to estimate Wooldridge #
prodestWRDG <- function(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL, seed = 123456,
                        tol = 1e-100, theta0 = NULL, cluster = NULL){
  set.seed(seed)
  Start = Sys.time() # start tracking time
  Y <- checkM(Y) # change all input to matrix
  fX <- checkM(fX)
  sX <- checkM(sX)
  pX <- checkM(pX)
  idvar <- checkM(idvar)
  timevar <- checkM(timevar)
  opt = 'optim' # this is going to be changed once the routine will work with DEoptim and solnp
  snum <- ncol(sX) # find the number of input variables
  fnum <- ncol(fX)
  if (!is.null(cX)) {cX <- checkM(cX); cnum <- ncol(cX)} else {cnum <- 0} # if is there any control, take it into account, else fix the number of controls to 0
  lag.fX = fX # generate fX lags
  for (i in 1:fnum) {
    lag.fX[, i] = lagPanel(fX[, i], idvar = idvar, timevar = timevar)
  }
  polyframe <- data.frame(sX,pX) # vars to be used in polynomial approximation
  regvars <- cbind(model.matrix( ~.^2-1, data = polyframe),sX^2,pX^2) # generate a polynomial of the desired level
  lagregvars <- regvars
  for (i in 1:dim(regvars)[2]) {
    lagregvars[, i] <- lagPanel(idvar = idvar, timevar = timevar, regvars[ ,i])
  }
  data <- model.frame(Y ~ fX + sX + lag.fX + regvars[, -1] + lagregvars + idvar + timevar) # data.frame of usable observations --> regvars
  if (is.null(theta0)) {
    theta0 <- runif((2 + fnum + snum + cnum + ncol(data$regvars) + ncol(data$lagregvars)),0,1)
  } else {
    if (length(theta0) != (fnum + snum)){ # WRONG NUMBER OF STARTING POINTS
      stop(paste0('theta0 length (', length(theta0), ') is inconsistent with the number of parameters (', fnum + snum, ')'), sep = '')
    } else { # CORRECT STARTING POINTS
      theta0 <- c(runif(1,0,1), theta0, runif(1 + cnum + ncol(data$regvars) + ncol(data$lagregvars)))
    }
  }
  vY <- data$Y
  X1 = cbind(1, data$fX, data$sX, data$regvars)
  X2 = cbind(1, data$fX, data$sX, 1, data$lagregvars)
  Z1 = cbind(1, data$fX, data$sX, data$regvars)
  Z2 = cbind(1, data$lag.fX, data$sX, data$lagregvars) # define the matrices of data + instruments

  betas.1st <- optim(theta0, gWRDG, method = "BFGS", Y = vY, X1 = X1, X2 =  X2, Z1 = Z1,
                 Z2 = Z2, numR = fnum + snum + cnum)
  W.star <- weightM(Y = vY, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, betas = betas.1st$par,
                    numR = (fnum + snum + cnum + 1)) # Compute the optimal weighting matrix
  betas <- optim(betas.1st$par, gWRDG, method = "BFGS", Y = vY, X1 = X1, X2 =  X2, Z1 = Z1,
                 Z2 = Z2, numR = fnum + snum + cnum, W = W.star)

  ### BOOTSTRAPPING THE STANDARD ERRORS --> THIS IS NOT THE LETTER OF WOOLDRIDGE (2009), TEMPORARY PART
  boot.indices <- block.boot.resample(data$idvar, R) # generates a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    boot.betas <- matrix(unlist(lapply(boot.indices, finalWRDG, data = data,
             theta0 = theta0, opt = opt, W = W.star)), ncol = fnum + snum + cnum, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  } else {
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix(unlist(
      parLapply(cl = cluster, boot.indices, finalWRDG, data = data,
                theta0 = theta0, opt = opt, W = W.star)), ncol = fnum + snum + cnum, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas, 2, sd, na.rm = TRUE) # calculate standard deviations

  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(cX)) {res.names <- c(res.names, colnames(cX, do.NULL = FALSE, prefix = 'cX'))} # add cX names to the results' names
  betapar <- betas$par[2 : (snum + fnum + cnum + 1)]
  betase <- boot.errors
  names(betapar) <- res.names # change results' names
  names(betase) <- res.names # change results' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'WRDG', boot.repetitions = NA, elapsed.time = elapsedTime, startingPoints = NA ,
                          opt = opt, seed = seed, opt.outcome = betas),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar),
             Estimates = list(pars = betapar, std.errors = betase))
  return(out)
}
# end of Wooldridge #

# function to prepare data for bootstrapping #
finalWRDG <- function(ind, data, idvar, timevar, theta0, opt, W){
  if (ind[1] != TRUE){
    idvar <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  } else {
    idvar <- data[ind, 'idvar', drop = FALSE]
  }
  data <- data[ind,] # change the index according to bootstrapped indices
  fnum <- ncol(data$fX)
  snum <- ncol(data$sX)
  Y <- data$Y
  X1 = cbind(1, data$fX, data$sX, data$regvars)
  X2 = cbind(1, data$fX, data$sX, 1, data$lagregvars)
  Z1 = cbind(1, data$fX, data$sX, data$regvars)
  Z2 = cbind(1, data$lag.fX, data$sX, data$lagregvars) # define the data and the matrices
  if (is.null(dim(data$cX))) {
    cnum <- 0
  } else {
    cX <- dim(data$cX)[2]
  }
  if (opt == 'optim'){
    try.out <- try(optim(theta0, gWRDG, method = "BFGS", Y = Y, X1 = X1, X2 = X2, Z1 = Z1,
                         Z2 = Z2, numR = fnum + snum + cnum, W = W), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      boot.beta <- try.out$par
    } else {
      boot.beta <- matrix(NA, (snum + fnum + cnum + 1), 1)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'solnp'){
    try.out <- try(solnp(theta0, gWRDG, Y = Y, X1 = X1, X2 =  X2, Z1 = Z1,
                        Z2 = Z2, numR = fnum + snum + cnum, W = W,
                        control = list(trace = FALSE)), silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      boot.beta <- try.out$pars
    } else {
      boot.beta <- matrix(NA, (snum + fnum + cnum + 1), 1)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'DEoptim'){
    try.out <- try(DEoptim(fn = gWRDG, lower = c(theta0), upper = rep.int(1,length(theta0)),
                          Y = Y, X1 = X1, X2 =  X2, Z1 = Z1, Z2 = Z2, numR = fnum + snum + cnum,
                          W = W, control = DEoptim.control(trace = FALSE)),
                  silent = TRUE)
    if (!inherits(try.out, "try-error")) {
      boot.beta <- try.out$optim$bestmem
    } else {
      boot.beta <- matrix(NA, (snum + fnum + cnum + 1), 1)
    } # error handling: if the optimization fails return missing values
  }
  return(boot.beta[2 : (fnum + snum + cnum + 1)])
}
# end of boot function #

# function to run the GMM estimation for WRDG #
gWRDG <- function(theta, Y, X1, X2, Z1, Z2, numR, W = NULL){
  k1 <- ncol(X1)
  R1t <- Y - X1 %*% theta[1 : k1, drop = FALSE]
  M1 <- t(Z1) %*% R1t
  R2t <- Y - X2 %*% c(theta[1:(numR + 1)], theta[(k1 + 1) : length(theta), drop = FALSE])
  M2 <- t(Z2) %*% R2t
  M <- rbind(M1, M2)
  if (is.null(W)){
    Z <- as.matrix(bdiag(Z1, Z2))
    W <- solve((t(Z) %*% Z) * ncol(Z)) # Initial matrix: unadjusted and independent --> diagonal
  }
  crit <- t(M) %*% W %*% M
  return(crit)
}
# end of GMM WRDG #
