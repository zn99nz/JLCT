source('JLCT\\JLCA_code_UD_2.r')
source('JLCT\\set_init.r')

my_JLCA_initset <- function(data, n, njoint, nc, nclass, nitem, weight = rep(1, n), 
                            init = NA, rho_fix = FALSE, rho_init = NA, 
                            init_iter = 25, init_num = 30,print = FALSE){
  
  c1 <- length(weight) != n
  if(c1){return('Please check weights.')}else{
    
    if(all(is.na(unlist(init)))){
      # select best initial value among 10 sets of initial value.
      log_like_temp <- 0
      # print('selecting best initial value')
      
      for(i in 1:init_num){
        init_temp <- set_init(njoint, nc, nclass, nitem)
        
        # rho_fix
        if(rho_fix){init_temp$rho <- rho_init}
        
        JLCA_init <- my_JLCA_UD(data, n, njoint, nc, nclass, nitem, weight,
                                init_temp, rho_fix, init_iter, print = FALSE)
        
        # update to better initial value 
        if((JLCA_init$log_like > log_like_temp) | (log_like_temp == 0)){
          log_like_temp <- JLCA_init$log_like
          init <- init_temp
          set <- i
        }
      }
      
      JLCA_result <- my_JLCA_UD(data, n, njoint, nc, nclass, nitem, weight, init, rho_fix, print = print)
    }else{
      JLCA_result <- my_JLCA_UD(data, n, njoint, nc, nclass, nitem, weight, init, rho_fix, print = print)
    } # if_else : setting initial value
  } # if_else : check for weight
  
  return(JLCA_result)
} # function

