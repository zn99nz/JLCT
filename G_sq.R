sourceCpp('C:\\Users\\uccoo\\Desktop\\ÇÐ±³\\study\\JLCT\\G_sq.cpp')

G_sq <- function(data, params){
  gamma <- params$gamma
  eta <- params$eta
  rho <- params$rho
  
  njoint <- length(gamma)
  nc <- length(eta)
  nclass <- sapply(eta, ncol)
  nitem <- ncol(data[[1]])
  n <- nrow(data[[1]])
  
  ex_like <- G_sq_ex(njoint, nitem, n, nc, nclass, gamma, eta, rho, data, weight=rep(1,n))
  
  comb_data <- data[[1]]
  for(i in 2:length(data)){
    comb_data <- cbind(comb_data, data[[i]])
  }
  
  rsp_data <- apply(comb_data, 1, paste, collapse='')
  unq_rsp <- unique(apply(comb_data, 1, paste, collapse=''))
  names(rsp_data) <- names(unq_rsp) <- NULL
  
  obs_f <- ex_f <- rep(0, length(unq_rsp))
  
  for(i in 1:length(unq_rsp)){
    obs_f[i] <- sum(unq_rsp[i]==rsp_data)
    pat_ind <- which(unq_rsp[i]==rsp_data)[1]
    ex_f[i] <- ex_like[pat_ind]*n
  }
  
  gsq <- 2*sum(obs_f*log(obs_f/ex_f))
  return(gsq)
}




