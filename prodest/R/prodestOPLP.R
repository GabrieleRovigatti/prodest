############ OLLEY-PAKES and LEVINSOHN-PETRIN ###############

# function to estimate OP and LP model #
prodestOP <- function(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL, opt = 'optim',
                      theta0 = NULL, seed = 123456, cluster = NULL, tol = 1e-100){
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
  polyframe <- data.frame(cbind(sX,pX)) # vars to be used in polynomial approximation
  mod <- model.matrix( ~ .^2 - 1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe), rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(fX, cX, mod, sX^2, pX^2) # generate a polynomial of the desired level
  lag.sX = sX # generate sX lags
  for (i in 1:snum) {
    lag.sX[, i] = lagPanel(sX[, i], idvar = idvar, timevar = timevar)
  }
  if (!is.null(cX)) {  # generate the matrix of data
    data <- suppressWarnings(as.matrix(data.frame(state = sX, lag.sX = lag.sX, free = fX, cX = cX, Y = Y,
                                 idvar = idvar, timevar = timevar, regvars = regvars)))
  } else {
    data <- suppressWarnings(as.matrix(data.frame(state = sX, lag.sX = lag.sX, free = fX, Y = Y, idvar = idvar,
                                 timevar = timevar, regvars = regvars)))
  }
  betas <- finalOPLP(ind = TRUE, data = data, fnum = fnum, snum = snum, cnum = cnum, opt = opt,
                     theta0 = theta0, boot = FALSE, tol = tol) # generate the list of estimated betas
  boot.indices <- block.boot.resample(idvar, R) # generate a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    nCores = NULL
    boot.betas <- matrix(unlist(lapply(boot.indices, finalOPLP, data = data, fnum = fnum, snum = snum, cnum = cnum,
                                       opt = opt, theta0 = theta0, boot = TRUE, tol = tol)), ncol = fnum+snum+cnum, byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  } else { # set up the cluster: send the lag Panel and the data.table libraries to clusters
    nCores = length(cluster)
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix(unlist(parLapply(cl = cluster, boot.indices, finalOPLP, data = data,
                                          fnum = fnum, snum = snum, cnum = cnum,
                                          opt = opt, theta0 = theta0, boot = TRUE, tol = tol)), ncol = fnum+snum+cnum,
                         byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas,2,sd,na.rm=TRUE) # calculate standard deviations
  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(cX)) {res.names <- c(res.names,colnames(cX, do.NULL = FALSE, prefix = 'cX'))}
  names(betas$betas) <- res.names # change results' names
  names(betas$FSbetas)[2 : length(betas$FSbetas)] <- c(res.names, rep(" ", (length(betas$FSbetas) - length(res.names) - 1))) # change 1st stage results' names
  names(boot.errors) <- res.names # change std.errors' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'OP', FSbetas = betas$FSbetas, boot.repetitions = R, elapsed.time = elapsedTime,
                          theta0 = theta0 , opt = opt, seed = seed, opt.outcome = betas$opt.outcome, nCores = nCores),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar,
                         FSresiduals = betas$FSresiduals),
             Estimates = list(pars = betas$betas, std.errors = boot.errors))
  return(out)
}
# end of prodestOP #

# function to estimate OP and LP model #
prodestLP <- function(Y, fX, sX, pX, idvar, timevar, R = 20, cX = NULL, opt = 'optim',
                      theta0 = NULL, seed = 123456, cluster = NULL, tol = 1e-100){
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
  polyframe <- data.frame(cbind(sX,pX)) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(fX, cX, mod, sX^2, pX^2) # generate a polynomial of the desired level
  lag.sX = sX # generate sX lags
  for (i in 1:snum) {
    lag.sX[, i] = lagPanel(sX[, i], idvar = idvar, timevar = timevar)
  }
  if (!is.null(cX)) {  # generate the matrix of data
    data <- suppressWarnings(as.matrix(data.frame(state = sX, lag.sX = lag.sX, free = fX, cX = cX, Y = Y,
                                 idvar = idvar, timevar = timevar, regvars = regvars)))
  } else {
    data <- suppressWarnings(as.matrix(data.frame(state = sX, lag.sX = lag.sX, free = fX, Y = Y, idvar = idvar,
                                 timevar = timevar, regvars = regvars)))
  }
  betas <- finalOPLP(ind = TRUE, data = data, fnum = fnum, snum = snum, cnum = cnum, opt = opt,
                     theta0 = theta0, boot = FALSE, tol = tol) # generate the list of estimated betas
  boot.indices <- block.boot.resample(idvar, R) # generate a list: every element has different length (same IDs, different time occasions) and is a vector of new indices, whose rownames are the new IDs
  if (is.null(cluster)){
    nCores = NULL
    boot.betas <- matrix(NA, R, (fnum + snum + cnum))
    for (i in 1:R){
      boot.betas[i,] <- finalOPLP(ind = boot.indices[[i]], data = data, fnum = fnum, snum = snum, cnum = cnum,
                                  opt = opt, theta0 = theta0, boot = TRUE, tol = tol)
    }
  } else { # set up the cluster: send the lag Panel and the data.table libraries to clusters
    nCores = length(cluster)
    clusterEvalQ(cl = cluster, library(prodest))
    boot.betas <- matrix(unlist(parLapply(cl = cluster, boot.indices, finalOPLP, data = data,
                                          fnum = fnum, snum = snum, cnum = cnum,
                                          opt = opt, theta0 = theta0, boot = TRUE, tol = tol)), ncol = fnum + snum + cnum,
                         byrow = TRUE) # use the indices and pass them to the final function (reshape the data)
  }
  boot.errors <- apply(boot.betas, 2, sd, na.rm = TRUE) # calculate standard deviations
  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(cX)) {res.names <- c(res.names,colnames(cX, do.NULL = FALSE, prefix = 'cX'))}
  names(betas$betas) <- res.names # change results' names
  names(betas$FSbetas)[2 : length(betas$FSbetas)] <- c(res.names, rep(" ", (length(betas$FSbetas) - length(res.names) - 1))) # change 1st stage results' names
  names(boot.errors) <- res.names # change std.errors' names
  elapsedTime = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'LP', FSbetas = betas$FSbetas, boot.repetitions = R, elapsed.time = elapsedTime,
                          theta0 = theta0 , opt = opt, seed = seed, opt.outcome = betas$opt.outcome, nCores = nCores),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar,
                         FSresiduals = betas$FSresiduals),
             Estimates = list(pars = betas$betas, std.errors = boot.errors))
  return(out)
}
# end of prodestLP #

# function to estimate and to bootstrap OP / LP #
finalOPLP <- function(ind, data, fnum, snum, cnum, opt, theta0, boot, tol){
  if (sum(as.numeric(ind)) == length(ind)){ # if the ind variable is not always TRUE
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  data <- data[ind,] # change the index according to bootstrapped indices
  first.stage <- lm(data[,'Y'] ~ data[,grepl('regvars', colnames(data))], na.action = na.exclude)
  fX <- data[,grepl('free', colnames(data)), drop = FALSE]
  phi <- fitted(first.stage) # generate the fitted values of the first stage
  beta.free <-  coef(first.stage)[2:(1 + fnum)]
  if (cnum != 0) { # save betas of cX if included
    beta.control <- coef(first.stage)[(2 + fnum):(1 + fnum + cnum)]
  } else {
    beta.control <- NULL
  }
  if (is.null(theta0)) { # if not specified, starting points are the first stage + normal noise
    theta0 <- coef(first.stage)[(2 + fnum):(1 + fnum + snum)]  + rnorm((snum), 0, 0.01)
  }
  phi <- phi - (fX %*% beta.free) # "clean" the phi from the effects of free vars
  newtime <- data[,'timevar', drop = FALSE]
  rownames(phi) <- NULL
  rownames(newtime) <- NULL
  lag.phi <- lagPanel(value = phi, idvar = newid, timevar = newtime) #  # lag fitted values
  res <- data[,'Y', drop = FALSE] - (fX %*% beta.free) # "clean" the outcome from the effects of free vars
  state = data[, grepl('state', colnames(data)), drop = FALSE]
  lag.sX = data[, grepl('lag.sX', colnames(data)), drop = FALSE]
  tmp.data <- model.frame(state ~ lag.sX + phi + lag.phi + res)
  if (opt == 'optim'){
    try.state <- try(optim(theta0, gOPLP, method = "BFGS", mX = tmp.data$state, mlX = tmp.data$lag.sX,
                           vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi, vres = tmp.data$res, stol = tol), silent = TRUE)
    if (!inherits(try.state, "try-error")) {
      beta.state <- try.state$par
      opt.outcome <- try.state
    } else {
      beta.state <- matrix(NA,(snum),1)
      opt.outcome <- list(convergence = 99)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'DEoptim') {
    try.state <- try(DEoptim(fn = gOPLP, lower = c(theta0), upper = rep.int(1,length(theta0)),
                             mX = tmp.data$state, mlX = tmp.data$lag.sX,
                             vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi,
                             vres =  tmp.data$res, stol = tol, control = DEoptim.control(trace = FALSE)),
                     silent = TRUE)
    if (!inherits(try.state, "try-error")) {
      beta.state <- try.state$optim$bestmem
      opt.outcome <- try.state
    } else {
      beta.state <- matrix(NA,(snum), 1)
      opt.outcome <- list(convergence = 99)
    } # error handling: if the optimization fails return missing values
  } else if (opt == 'solnp') {
    try.state <- try(solnp(theta0, gOPLP, mX = tmp.data$state, mlX = tmp.data$lag.sX,
                           vphi = tmp.data$phi, vlag.phi = tmp.data$lag.phi,
                           vres =  tmp.data$res, stol = tol, control = list(trace = FALSE)), silent = TRUE)
    if (!inherits(try.state, "try-error")) {
      beta.state <- try.state$pars
      opt.outcome <- try.state
    } else {
      beta.state <- matrix(NA,(snum),1)
      opt.outcome <- list(convergence = 99)
    } # error handling: if the optimization fails return missing values
  }
  if (boot == FALSE){ # for the baseline estimation we want info on optimization, too
    return(list(betas = c(beta.free, beta.state, beta.control), opt.outcome = opt.outcome, FSbetas = coef(first.stage),
                FSresiduals = resid(first.stage)))
  } else{
    return(c(beta.free, beta.state, beta.control))
  }
}
# end of OP / LP final function #

# function to run the GMM estimation for OP and LP #
gOPLP <- function(vtheta, mX, mlX, vphi, vlag.phi, vres, stol){
  Omega <- vphi - mX %*% vtheta
  Omega_lag <- vlag.phi - mlX %*% vtheta
  Omega_lag_pol <- cbind(1, Omega_lag, Omega_lag^2, Omega_lag^3)
  g_b <- solve(crossprod(Omega_lag_pol), tol = stol) %*% t(Omega_lag_pol) %*% Omega
  XI <- vres - (mX %*% vtheta) - (Omega_lag_pol %*% g_b)
  crit <- crossprod(XI)
  return(crit)
}
# end of GMM OPLP #
