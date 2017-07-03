# function to demean a vector and take the variance #
withinvar <- function(inmat){
  devmat = inmat - mean(inmat) # demean the vector
  return(var(c(devmat)))
}
# end of within variance function #

## Simulate a panel dataset  ##
panelSim <- function(N = 1000, T = 100, alphaL = .6, alphaK = .4, DGP = 1, rho = .7, sigeps = .1,
                     sigomg = .3, rholnw = .3, seed = 123456){
  set.seed(seed)
  if (DGP != 1 & DGP != 2 & DGP != 3) stop('DGP must be 1, 2 or 3')
  if (T < 15) stop('T must be > 15 to yield a usable number of observations - i.e., 2 per panel.')
  siglnw.tot <- c(.1,0,.1) # parameters describing the different DGP
  timeb.tot <- c(.5,0,.5)
  sigoptl.tot <- c(0,.37,.37)
  measerrvec <- c(0,sqrt(.1),sqrt(.2),sqrt(.5)) # measurement error of the intermediate input
  starttime <- round(T*.9)
  numT <- T - starttime #this is the number of times to be used
  alpha0 <- 0 # starting points
  epsilon <- matrix(rnorm(N*T,0,sigeps),nrow=N,ncol=T) # matrix of random normal shocks
  siglnw <- siglnw.tot[DGP]
  timeb <- timeb.tot[DGP]
  sigoptl <- sigoptl.tot[DGP]
  # The next five lines define parameters necessary to subdivide the omega process into omega(t-1) -> omega(t-b) -> omega(t).  Note that this depends on the timeb (b) parameter, defined at top of program
  rhofirst <- rho^(1-timeb)
  rhosecond <- rho^timeb
  sigxi <- sqrt((1-rho^2)*sigomg^2)  # SD of innovation in omega - set such that variance of omega(t) is constant across time
  sigxifirst <- sqrt((1-rhofirst^2)*sigomg^2)
  sigxisecond <- sqrt((1-rhosecond^2)*sigomg^2)
  sigxilnw <- sqrt((1-rholnw^2)*siglnw^2) # SD of innovation in ln(wage) - set such that variance of ln(wage) is constant across time
  omgdata0 <- matrix(rnorm(N,0,sigomg),nrow=N,ncol=1)
  lnwdata0 <- matrix(rnorm(N,0,siglnw),nrow=N,ncol=1)
  # Period 1-b value of omega and Period 1 values of omega and ln(wage)
  omgdataminusb <- matrix(NA,N,T)
  omgdata <- matrix(NA,N,T)
  lnwdata <- matrix(NA,N,T)
  lnkdata <- matrix(NA,N,T)
  lnldata <- matrix(NA,N,T)
  lnpdata <- matrix(0,N,T)
  epsdata <- matrix(rnorm(N*T,0,sigeps),N,T)
  omgdataminusb[,1] = rhofirst*omgdata0 + matrix(rnorm(N,0,sigxifirst),nrow=N,ncol=1)
  omgdata[,1] = rhosecond*omgdataminusb[,1] + matrix(rnorm(N,0,sigxisecond),nrow=N,ncol=1)
  lnwdata[,1] = rholnw*lnwdata0 + matrix(rnorm(N,0,sigxilnw),nrow=N,ncol=1)
  # Simulate values of omega and ln(wage) for rest of time periods
  for (i in 1:N) {
    for (t in 2:T) {
      omgdataminusb[i,t] = rhofirst*omgdata[i,t-1] + rnorm(1,0,sigxifirst)
      omgdata[i,t] = rhosecond*omgdataminusb[i,t] + rnorm(1,0,sigxisecond)
      lnwdata[i,t] = rholnw*lnwdata[i,t-1] + rnorm(1,0,sigxilnw)
    }
  }
  lnkdata[,1] = matrix(-100,N,1)
  disc = 0.95 # Discount rate for dynamic programming problem
  delta = 0.2 # Depreciation rate for capital
  sigb = 0.6 # measures variation in capital adjustment cost across firms
  oneoverbiadj = exp(rnorm(N,0,sigb))
  # Various components of long expression for optimal investment choice (end of Appendix Section 7.3)
  squarebracketterm = (alphaL^(alphaL/(1-alphaL)))*exp(0.5*alphaL^2*sigoptl^2) - (alphaL^(1/(1-alphaL)))*exp(0.5*sigoptl^2)
  const1 = disc * (alphaK/(1-alphaL)) * (exp(alpha0)^(1/(1-alphaL))) * squarebracketterm
  vec1 = (disc*(1-delta))^seq(100)
  vec2 = cumsum( rholnw^(2*seq(100)) ) # Cumulative sum done through an upper triangular matrix --> check
  vec3 = rbind(sigxi^2 * 0,cumsum(rho^(2*(seq(100)-1))) )
  expterm3 = exp( 0.5 * ((-alphaL)/(1-alphaL))^2 * sigxilnw^2 * vec2 );
  expterm4 = exp( 0.5 * (1/(1-alphaL))^2 * rhosecond^2 * ((sigxifirst^2)*rho^(2*seq(100)) + vec3) )
  expterm5 = exp((1/(1-alphaL))*(1/2)*sigxisecond^2)
  investmat = matrix(NA,N,T)
  for (i in 1:N) {
    for (t in 1:T) {
      expterm1 = exp( (1/(1-alphaL))*omgdata[i,t]*(rho^seq(100)) ) # first term in exponent in second line
      expterm2 = exp( ((-alphaL)/(1-alphaL))*lnwdata[i,t]*(rholnw^seq(100)) ) # second term in exponent in second line
      investmat[i,t] = oneoverbiadj[i]*const1*expterm5*sum(vec1*expterm1*expterm2*expterm3*expterm4)  # optimal investment
      if (t >= 2) {
        lnkdata[i,t] = log( (1-delta)*exp(lnkdata[i,t-1]) + (1-0*runif(1))*investmat[i,t-1] )
      }
    }
  }
  # Now generate levels of labor input for all firms and time periods - note: this choice depends on omega(t-b) since labor chosen at t-b
  for (i in 1:N) {
    for (t in 1:T) {
      lnldata[i,t] = ((sigxisecond^2)/2 + log(alphaL) + alpha0 + rhosecond*omgdataminusb[i,t] - lnwdata[i,t] + lnpdata[i,t] + (alphaK)*lnkdata[i,t])/(1-alphaL)
    }
  }
  # Now add potential optimization error to ln(labor) - note: sigoptl defined at top of program
  truelnldata = lnldata
  lnldata = lnldata + matrix(rnorm(N*T,0,sigoptl),N,T)
  # Generate levels of Output and Materials
  lnydata = alphaL*lnldata + alphaK*lnkdata + omgdata + epsdata
  lnmdata = alphaL*truelnldata + alphaK*lnkdata + omgdata
  truelnmdata = lnmdata
  withinvarlnm = withinvar(lnmdata[,(starttime+1):T]);
  # Note: Material input based on labor _without_ optimization error - this is necessary for LP to produce consistent estimates under DGP2 (without measurement error in materials) - see discussion in text. One interpretation of this is that econometrician observed planned (or ordered) materials (rather than materials actually used) */
  # Now add potential measurement error in material input - note that as described in the paper and table the amount is proportional to the within variance of ln(materials) */
  lnmdata2 = lnmdata + measerrvec[2]*matrix(rnorm(N*T,0,sqrt(withinvarlnm)),N,T)
  lnmdata3 = lnmdata + measerrvec[3]*matrix(rnorm(N*T,0,sqrt(withinvarlnm)),N,T)
  lnmdata4 = lnmdata + measerrvec[4]*matrix(rnorm(N*T,0,sqrt(withinvarlnm)),N,T)
  info <- list(lnydata,lnkdata,lnldata,lnmdata,lnmdata2,lnmdata3, lnmdata4)
  data <- matrix(unlist(lapply(info, function(x) { return(c(t((x[,(starttime+1):T])))) } )), ncol = 7)
  data <- data.frame(idvar = rep(1:N, each = numT), timevar = rep(1:numT,N), Y = data[,1], sX = data[,2], fX = data[,3],
                     pX1 = data[,4], pX2 = data[,5], pX3 = data[,6], pX4 = data[,7])
  return(data)
}



