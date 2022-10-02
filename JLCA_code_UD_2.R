library(Rcpp)
#############################################################################
#############################################################################

# equal number of item for each C latent variable
# All item is binary

sourceCpp('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\posterior_UD.cpp')

my_JLCA_UD <- function(data, n, njoint, nc, nclass, nitem, weight = rep(1, n), 
                       init, rho_fix = FALSE, maxiter = NA, print = TRUE){
  # njoint : # of joint class
  # nclass : # of latent class for each latent class variable
  # nc : # of latent class variable
  # weight : weight for JLCT (for JLCA, weight = 1/n)
  
  # data : list ((n*nitem)*nc)
  
  # initial value
  gamma <- init$gamma
  eta <- init$eta
  rho <- init$rho
  
  # posterior with initial value
  posterior <- post_UD(njoint, nitem, n, nc, nclass,
                       gamma, eta, rho, data, weight)
  
  post_jc <- posterior[[1]]
  post_j <- posterior[[2]] # prob(joint class) : njoint*n
  post_c <- list()
  for(c in 1:nc){
    mat <- sapply(sapply(post_jc, '[', c), apply, 1, sum)
    post_c[[c]] <- mat
  }
  # prob(latent class) : list ((nclass*n)*nc)
  
  ########### start loop ###########
  iter <- 0; diff <- 1
  while(diff > 1e-4){
    iter <- iter+1
    
    # parameter update
    gamma_old <- gamma
    eta_old <- eta
    rho_old <- rho
    
    ##################################################################################3
    
    # 1) parameter update
    # 1-1) gamma
    w.mat <- matrix(weight, njoint, n, byrow = T)
    gamma <- as.vector(post_j%*%weight)/sum(weight)
    
    
    # 1-2) eta
    eta <- list()
    din <- as.vector(post_j%*%weight)
    for(c in 1:nc){
      din_mat <- matrix(din, njoint, nclass[c])
      post_class <- sapply(post_jc, '[', c)
      eta_temp <- matrix(0, nclass[c], njoint)
      for(i in 1:n){
        eta_temp <- eta_temp+(post_class[[i]]*weight[i])
      }
      eta[[c]] <- t(eta_temp)/din_mat
    }
    
    
    # 1-3) rho
    rho <- list()
    if(rho_fix){rho <- rho_old}else{
      for(c in 1:nc){
        data_c <- ifelse(data[[c]] == 1, 1, 0)
        mis_num <- apply(is.na(data_c), 2, sum)
        
        w.mat <- matrix(weight, nclass[c], n, byrow = T)
        post_c1 <- post_c[[c]]*w.mat
        rho_temp <- matrix(0, nclass[c], nitem)
        
        for(cl in 1:nclass[c]){
          data_mis <- data_c
          if(any(mis_num != 0)){
            data_mis[is.na(data_c)] <- rep(rho_old[[c]][cl,], mis_num)
          }
          rho_temp[cl,] <- t(post_c1[cl,])%*%data_mis
        }
        din_mat <- matrix(apply(post_c1, 1, sum), nclass[c], nitem)
        rho[[c]] <- rho_temp/din_mat
      }
    }
    
    
    
    # decide whether stop iteration
    diff <- max(abs(c(as.vector(gamma-gamma_old), 
                      as.vector(unlist(eta)-unlist(eta_old)), 
                      as.vector(unlist(rho)-unlist(rho_old)))))
    
    
    ###########################################################################################################
    
    # 2) posterior
    posterior <- post_UD(njoint, nitem, n, nc, nclass,
                         gamma, eta, rho, data, weight = rep(1,n))
    
    post_jc <- posterior[[1]]
    post_j <- posterior[[2]] # prob(joint class) : njoint*n
    post_c <- list()
    for(c in 1:nc){
      mat <- sapply(sapply(post_jc, '[', c), apply, 1, sum)
      post_c[[c]] <- mat
    }
    log_like <- posterior[[3]]
    # prob(latent class) : list ((nclass*n)*nc)
    
    if(print){
      print(paste('iter :', iter, '/ diff =', round(diff, 10), '/ log-like =', round(log_like, 3))) 
    }
    
    # diff 문제
    if(is.nan(diff)){
      gamma <- gamma_old; eta <- eta_old; rho <- rho_old
      break
    }
    
    # maxiter가 있는 경우
    if(!is.na(maxiter)){
      if(iter == maxiter) break
    }
    
  } # while
  
  
  # for result
  ## sorting
  sort_ind <- order(gamma, decreasing = TRUE)
  gamma <- gamma[sort_ind]
  eta <- lapply(eta, '[', sort_ind,)
  post_j <- post_j[sort_ind,]
  
  ## gamma
  names(gamma) <- paste0('joint ', 1:njoint)
  
  ## eta
  names(eta) <- paste0('latent var ', 1:nc)
  if(njoint > 1){
    for(c in 1:nc){
      dimnames(eta[[c]]) <- list(paste0('joint ', 1:njoint), paste0('class ', 1:nclass[c]))
    }
  }
  
  ## rho
  names(rho) <- paste0('latent var ', 1:nc)
  for(c in 1:nc){
    dimnames(rho[[c]]) <- list(paste0('class ', 1:nclass[c]), paste0('item ', 1:nitem))
  }
  
  ## result
  params <- list(gamma = gamma, eta = eta, rho = rho)
  g <- njoint-1; e <- njoint*(sum(nclass-1)); r <- sum(nitem*nclass)
  P <- g+e+r
  BIC <- -2*log_like + P*log(n)
  result <- list(params = params, iteration = iter, 
                 log_like = log_like, BIC = BIC, post = post_j)
  
  return(result)
  
} # function




AIC <- function(jlca, n){
  loglike <- jlca$log_like
  P <- (jlca$BIC+(2*loglike))/log(n)
  
  aic <- -2*loglike + 2*P
  return(aic)
}