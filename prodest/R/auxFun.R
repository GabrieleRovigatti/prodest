# function to check and change input to matrix #
checkM <- function(input){ # , inputname = NA
  if (!is.matrix(input)) {
    out <- as.matrix(input)
  } else{
    out <- input
  }
  colnames(out) <- NULL
  return(out)
}
# end check matrix fun #

# function to lag variables within a panel #
lagPanel <- function( idvar, timevar, value, suffix = '' ){
  df <- data.frame( idvar, timevar, value)
  last.time <- df %>% filter(!is.na(timevar)) %>%
    mutate(timevar = timevar + 1, lagged_value = value, value = NULL)
  out <- as.matrix(df %>% left_join(last.time, by = c("idvar", "timevar")))[,4]
  colnames( out ) <- NULL
  return( out )
}
# end of lag panel #

# boot resampling on IDs: bootstrapping on individuals #
block.boot.resample <- function( idvar, R ){
  unique.ids <- unique(idvar) # find the unique values of panels in order to reshape the data
  panel.time.indices <- apply(unique.ids, 1, function(x) {return(list(which(idvar == x)))}) # find the time indices for each panel
  seq.indices <- 1:length(unique.ids) # the panel.time.indices list is indexed with sequential numbers: we mimic it
  boot.panel.id <- replicate(R, sample(seq.indices, replace = TRUE)) # generate a matrix of new IDs - R times
  new.indices <- list() # generate the matrix of the new indices
  ind <- 1:length(unique.ids)
  for (r in 1:R){ # for each boot rep we generate a vector of indices with rownames equal to a new - and fake - ID
    new.indices[[r]] <- cbind(unlist(mapply(function(x,y) {
      names(panel.time.indices[[x]][[1]]) <- rep(y,length(panel.time.indices[[x]][[1]]))
      return(list(panel.time.indices[[x]][[1]]))
    }, boot.panel.id[,r], ind))) # return a fake ID (sequential number) as row name and the index referring to the true ID
  }
  return(new.indices)
}
# end of block bootstrap function #

# function to compute the weighting matrix #
weightM <- function(Y, X1, X2, Z1, Z2, betas, numR, SE = FALSE){
  k1 <- ncol(X1)
  N <- nrow(X1)
  R1t <- Y - X1 %*% betas[1 : k1, drop = FALSE]
  R2t <- Y - X2 %*% c(betas[1 : numR], betas[(k1 + 1) : length(betas), drop = FALSE]) # (fnum + snum + cnum + 1)
  u <- c(R1t, R2t) # alternative, still working
  Z <- as.matrix( bdiag(Z1, Z2)) # drop the collinear constant
  sigma.rs <- (t(u) %*% u)
  S <- sigma.rs[1] * ( ( t(Z) %*% Z) ) # /N
  if (SE == TRUE){
    dX <- rbind( cbind( X1, matrix(0, N, (ncol(X2) - numR) ) ),
                 cbind( X2[, 1 : numR], matrix(0 , N, (ncol(X1) - numR) ), X2[,(numR + 1) : ncol(X2)]) ) # generate a "quasi-block" matrix with common columns NON-BLOCK
    var.beta <- (1/N) * solve( ( t(dX) %*% Z ) %*% solve(S) %*% (t(Z) %*% dX) ) # compute varCovar matrix
    st.errors <- sqrt(diag(var.beta))
    return(st.errors)
  } else {
    W = solve(S)
    return(W)
  }
}
# end of weighting matrix function #

# function to print lateX table of results #
printProd <- function(mods, modnames = NULL, parnames = NULL, outfile = NULL, ptime = FALSE, nboot = FALSE){
  if (!is.null(outfile)) (sink(outfile)) # write on a text file
  numMods <- length(mods)
  numPars <- length(mods[[1]]@Estimates$pars)
  cat(paste('\\begin{tabular}{', paste(rep('c',(numMods*2+1)), collapse = ''),'}',
            '\\hline\\hline', sep = '')) # print tabular header
  nm <- '\n'
  obs <- '\nN'
  time <- '\nTime'
  boot <- '\nBootRep'
  for (m in 1:numMods){ # generate first and last row: names (methods or user-supplied) and observations
    if (is.null(modnames)){
      nm <- paste(nm, mods[[m]]@Model$method, sep = ' & & ')
    }else{
      nm <- paste(nm, modnames[m], sep = ' & & ')
    }
    obs <- paste(obs, length(mods[[m]]@Data$Y), sep = ' & & ')
    time <- paste(time, round(mods[[m]]@Model$elapsed.time[[1]], digits = 2), sep = ' & & ')
    boot <- paste(boot, mods[[m]]@Model$boot.repetitions, sep = ' & & ')
  }
  nm <- paste(nm, '\\\\\\hline')
  obs <- paste(obs, '\\\\\\hline\\hline')
  cat(nm)
  for (p in 1:numPars){ # generate the table body row by row: names (vars or user-supplied),
    if (is.null(parnames)){
      betas <- paste('\n', names(mods[[1]]@Estimates$pars)[p])
    }
    else{
      betas <- paste('\n', parnames[p])
    }
    sigmas <- '\n'
    blank <- '\n'
    for (m in 1:numMods){
      betas <- paste(betas, round(mods[[m]]@Estimates$pars[p],digits = 3), sep = ' & & ')
      sigma <- paste('(', round(mods[[m]]@Estimates$std.errors[p],digits = 3), ')', sep = '')
      sigmas <- paste(sigmas, sigma , sep = ' & & ')
      blank <- paste(blank, ' & ',  sep = '')
    }
    betas <- paste(betas, '\\\\')
    sigmas <- paste(sigmas, '\\\\')
    blank <- paste(blank, '\\\\')
    cat(betas)
    cat(sigmas)
    cat(blank)
  }
  cat(blank)
  if (ptime == TRUE) (cat(paste(time, '\\\\', sep = '')))
  if (nboot == TRUE) (cat(paste(boot, '\\\\', sep = '')))
  cat(obs)
  cat('\n\\end{tabular}')
  if (!is.null(outfile)) (sink())
}
# end of latex print table #
