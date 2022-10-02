library(gtools) # for rdirichlet()

set_init <- function(njoint, nc, nclass, nitem){
  # gamma : (njoint)
  gamma <- as.vector(rdirichlet(1, rep(1, njoint)))
  
  # eta : list ((njoint * nclass) * nc)
  eta <- list()
  for(c in 1:nc){
    eta[[c]] <- rdirichlet(njoint, rep(1, nclass[c]))
  }
  
  # rho : list ((nclass*nitem) * nc)
  rho <- list()
  for(c in 1:nc){
    rho[[c]] <- matrix(runif(nclass[c]*nitem), nclass[c], nitem)
  }
  
  init_value <- list()
  init_value$gamma <- gamma 
  init_value$eta <- eta
  init_value$rho <- rho
  
  return(init_value)
}