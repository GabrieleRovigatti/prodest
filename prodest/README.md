# prodest
Production Function Estimation in R

`prodest` ([Rovigatti, 2017](https://CRAN.R-project.org/package=prodest)) implements most of the methods for production function estimation. The prodest package provides functions to simulate panel data, estimate the productivity, state and free variables parameters. Full description of the models, their issues and characteristics is available in [Mollisi et al., 2017](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2916753).

The latest stable version of `prodest` is available on CRAN at 'https://cran.r-project.org/package=prodest'.

The latest development version of `prodest` is available at 'https://github.com/GabrieleRovigatti/prodest/tree/master/prodest'.

Please cite `prodest` in publications:

Rovigatti, G. (2017).  
_prodest: Production Function Estimation in R_.  
R package.  

Rovigatti, G. (2017).  
Production Function Estimation in R: The prodest Package.  
Working paper.      
 
Installation: To install the latest stable version of `prodest` in R type `install.packages("prodest")`.

Example: In R, type 

  require(prodest)
  
  data(chilean) # Chilean data on production.
  
  #we fit a model with two free (skilled and unskilled), one state (capital) and one proxy variable (electricity) with two different optimizers
  
  LP.fit <- prodestLP(chilean$Y, fX = cbind(chilean$fX1, chilean$fX2), chilean$sX,
                        chilean$pX, chilean$idvar, chilean$timevar, seed = 154673)
                        
  LP.fit.solnp <- prodestLP(chilean$Y, fX = cbind(chilean$fX1, chilean$fX2), chilean$sX,
                        chilean$pX, chilean$idvar, chilean$timevar, opt = 'solnp')
                        
  #show results
  
  summary(LP.fit)
  
  summary(LP.fit.solnp)

  #show results in .tex tabular format
  
  printProd(list(LP.fit, LP.fit.solnp))
