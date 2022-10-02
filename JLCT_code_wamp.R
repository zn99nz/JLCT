source('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\JLCA_code_initset_2.r')
source('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\JLCT_full_split.r')

JLCT <- function(data, n, nc, nclass, nitem, rho_fix = FALSE, rho_init = NA){
  
  print('[level 0]')
  print('split? JC_1 & JC_2')
  # for first split
  case <- combn(nc,2); best_wamp <- 0
  for(i in 1:ncol(case)){
    print(paste0(' * JLCA', i))
    t_case <- case[,i]
    t_data <- data[t_case]
    
    rho_init_t <- rho_init[t_case]
    t_JLCA <- my_JLCA_initset(t_data, n, 2, 2, nclass[t_case], 
                              nitem, weight = rep(1, n), 
                              rho_fix = rho_fix, rho_init = rho_init_t)
    t_wamp <- WAMP(t_JLCA$post, rep(1,n))
    
    if(best_wamp < t_wamp){
      best_wamp <- t_wamp
      best_case <- t_case
      best_JLCA <- t_JLCA
    }
  }
  print('yes')
  print('-------')
  
  # save the split result
  split_result <- list()
  
  
  for(i in 1:2){
    # parameter
    gamma <- best_JLCA$params$gamma[i]
    eta <- lapply(best_JLCA$params$eta, '[', i,)
    params <- list(gamma = gamma, eta = eta, rho = best_JLCA$params$rho)
    
    # weight
    weight <- best_JLCA$post[i,]
    
    # result
    split_result[[i]] <- list(var = best_case, params = params, weight = weight)
  }
  
  
  names(split_result) <- c('JC_1', 'JC_2')
  result <- list('level_1' = split_result)
  level <- 0
  ########################################################
  node_split <- TRUE
  while(length(result) > level){
    level <- level+1
    
    print(paste0('[level ', level, ']'))
    result <- JLCT_full_split(data, n, nclass, nitem, result, level, rho_fix, rho_init)
    
  }
  
  return(result)
} # function


