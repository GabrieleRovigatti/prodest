############ WOOLDRIDGE ###############

# function to estimate Wooldridge #
prodestWRDG_GMM <- function(Y, fX, sX, pX, idvar, timevar, cX = NULL, seed = 123456, tol = 1e-100){
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
  data <- model.frame(Y ~ fX + sX + lag.fX + regvars + lagregvars + idvar + timevar) # data.frame of usable observations --> regvars

  Y <- data$Y; dY = c(data$Y, data$Y); X1 = cbind(1, data$fX, data$regvars)
  X2 = cbind(1, data$fX, data$sX, 1, data$lagregvars)
  Z1 = cbind(1, data$fX, data$regvars)
  Z2 = cbind(1, data$lag.fX, data$sX, data$lagregvars) # define the data and the matrices
  numR = 1 + fnum + snum + cnum # number of restricted columns
  numU1 <- ncol(X1) - numR # number of unrestricted cols of X1
  numU2 <- ncol(X2) - numR # number of unrestricted cols of X2
  N <- nrow(X1)

  dX <- rbind( cbind( X1, matrix(0, N, numU2 ) ),
               cbind( X2[, 1 : numR], matrix(0 , N, numU1), X2[,(numR + 1) : ncol(X2)]) ) # generate a "quasi-block" matrix with common columns NON-BLOCK
  Z <- as.matrix(bdiag(Z1, Z2))
  W <- solve((t(Z) %*% Z)) * diag(ncol(Z)) # unadjusted, independent

  betas.1st <- solve(t(dX) %*% Z %*% W %*% t(Z) %*% dX, tol = tol) %*%
                  t(dX) %*% Z %*% W %*% t(Z) %*% dY # 1st stage GMM parameters
  W.star <- weightM(Y = Y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
                    betas = betas.1st, numR = (fnum + snum + cnum + 1)) # compute optimal weighting matrix
  betas.2nd <- solve(t(dX) %*% Z %*% W.star %*% t(Z) %*% dX, tol = tol) %*%
                  t(dX) %*% Z %*% W.star %*% t(Z) %*% dY # 2nd step
  st.errors <- weightM(Y = Y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
                       betas = betas.2nd, numR = (fnum + snum + cnum + 1), SE = TRUE) # compute standard errors

  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results
  if (!is.null(cX)) {res.names <- c(res.names, colnames(cX, do.NULL = FALSE, prefix = 'cX'))} # add cX names to the results' names
  betapar <- betas.2nd[2: (snum + fnum + cnum + 1)]
  betase <- st.errors[2: (snum + fnum + cnum + 1)]
  names(betapar) <- res.names # change results' names
  names(betase) <- res.names # change results' names
  elapsed.time = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'WRDG', boot.repetitions = NA, elapsed.time = elapsed.time, theta0 = NA,
                          opt = NA, seed = seed, opt.outcome = NULL),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar),
             Estimates = list(pars = betapar, std.errors = betase))
  return(out)
}
# end of Wooldridge #
