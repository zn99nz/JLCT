# full split
# 각 node에서 비교해야하는 JLC model을 끝까지 적합 후 best JLC model을 결정 및 split


#entropy <- function(post, weight){
#  n <- nrow(post); C <- ncol(post)
#  s_post <- sum(-weight*post*log(post))
#  ent <- 1-(s_post/(sum(weight)*log(C)))
#  return(ent)
#}


WAMP <- function(post, weight){
  max_post <- apply(post, 1, max)
  wamp <- sum(weight*post)/sum(weight)
  return(wamp)
}

JLCT_full_split <- function(data, n, nclass, nitem, result, level, rho_fix = FALSE, rho_init = NA){
  # present level (= parent nodes)
  p_result <- result[[level]]
  
  # for saving split result (= child nodes split from each parent node)
  c_result <- list()
  
  # for each node
  # parent node 갯수
  np_node <- length(p_result)
  for(nd in 1:np_node){
    # 현재 child node 갯수
    nc_node <- length(c_result)
    
    p_node <- p_result[[nd]]; p_name <- names(p_result)[nd]
    print(paste('split?', paste0(p_name, c(1,2))[1], '&', paste0(p_name, c(1,2))[2]))
    p_var <- p_node$var; p_data <- data[p_var]
    weight <- p_node$weight
    
    BIC_list <- NULL
    WAMP_list <- NULL
    
    ## 2-joint class JLCA with existing variables
    print(' * JLCA_0')
    
    rho_init1 <- rho_init[p_var]
    
    JLCA_1 <- my_JLCA_initset(p_data, n, njoint=2, nc=length(p_var), 
                              nclass=nclass[p_var], nitem, weight, 
                              rho_fix = rho_fix, rho_init = rho_init1)
    
    best_JLCA <- JLCA_1
    best_wamp <- WAMP(best_JLCA$post, weight)
    best_BIC <- best_JLCA$BIC
    WAMP_list <- c(WAMP_list, best_wamp)
    BIC_list <- c(BIC_list, best_BIC)
    
    
    best_case <- 0
    
    ## JLCA adding a new variable
    ## 고려하고 있는 node가 모든 class var를 사용한 경우는 pass
    jlca <- 0
    if(length(p_var) < length(data)){
      for(v in 1:length(data)){
        if(!(v %in% p_var)){
          jlca <- jlca+1
          print(paste0(' * JLCA+', jlca))
          t_var <- sort(c(p_var, v))
          t_data <- data[t_var]
          
          rho_init_t <- rho_init[t_var]
          t_JLCA <- my_JLCA_initset(t_data, n, njoint=2, nc=length(t_var), 
                                    nclass=nclass[t_var], nitem, weight,
                                    rho_fix = rho_fix, rho_init = rho_init_t)
          
          t_wamp <- WAMP(t_JLCA$post, weight)
          WAMP_list <- c(WAMP_list, t_wamp)
          BIC_list <- c(BIC_list, t_JLCA$BIC)
          
          if(jlca == 1){
            t_best_wamp <- WAMP(t_JLCA$post, weight)
            t_best_JLCA <- t_JLCA
            t_best_case <- v
          }else{
            if(t_wamp > t_best_wamp){
              t_best_wamp <- t_wamp
              t_best_JLCA <- t_JLCA
              t_best_case <- v
            }
          } # JLCA 결과 업데이트
        } # if !(v %in% p_var)
      } # for v
      
      if(best_wamp < t_best_wamp){
        best_JLCA <- t_best_JLCA
        best_case <- t_best_case
      }
      
      #t_best_ent <- entropy(t(t_best_JLCA$post), weight)
      
      #if(best_ent < t_best_ent){
      #  best_JLCA <- t_best_JLCA
      #  best_case <- t_best_case
      #}
    }
    
    
    ## 1-joint class JLCA with existing variables
    if(best_case == 0){
      print(' * JLCA-1')
      
      JLCA_2 <- my_JLCA_initset(p_data, n, njoint=1, nc=length(p_var), 
                                nclass=nclass[p_var], nitem, weight, 
                                rho_fix = rho_fix, rho_init = rho_init1)
      BIC_list <- c(BIC_list, JLCA_2$BIC)
      
      if(JLCA_2$BIC < best_BIC){
        best_JLCA <- JLCA_2
        best_case <- -1
        print('no')
      }else{print('yes')}
    }else{print('yes')}
    
    
    ######################## determine the best JLC model
    
    ## best model -> split or not
    p_result[[nd]]$BIC <- BIC_list
    p_result[[nd]]$WAMP <- WAMP_list
    
    if(best_case == -1){
      # not split
      p_result[[nd]]$split <- FALSE

    }else{
      if(best_case == 0){
        # split using existing variables
        p_result[[nd]]$split <- TRUE
        
        for(i in 1:2){
          # parameter
          gamma <- best_JLCA$params$gamma[i]
          eta <- lapply(best_JLCA$params$eta, '[', i,)
          params <- list(gamma = gamma, eta = eta, rho = best_JLCA$params$rho)
          
          # weight
          c_weight <- best_JLCA$post[i,]*weight
          
          # result
          c_result[[nc_node+i]] <- list(var = p_var, params = params, weight = c_weight)
        }
        names(c_result)[nc_node+(1:2)] <- paste0(p_name, 1:2)
        
      }else{
        # split adding a new variable
        p_result[[nd]]$split <- TRUE
        
        for(i in 1:2){
          # parameter
          gamma <- best_JLCA$params$gamma[i]
          eta <- lapply(best_JLCA$params$eta, '[', i,)
          params <- list(gamma = gamma, eta = eta, rho = best_JLCA$params$rho)
          
          # weight
          c_weight <- best_JLCA$post[i,]*weight
          
          # result
          c_result[[nc_node+i]] <- list(var = sort(c(p_var, best_case)),
                                        params = params, weight = c_weight)
        }
        
        names(c_result)[nc_node+(1:2)] <- paste0(p_name, 1:2)
      }
    }
  } # for nd
  
  result[[level]] <- p_result # each parent node의 split 정보 update
  
  # split된 child node 정보 update
  if(length(c_result) > 0){
    result[[level+1]] <- c_result
    names(result)[level+1] <- paste('level_', level+1, sep = '')
  }
  
  print('-------')
  
  return(result)
} # function


