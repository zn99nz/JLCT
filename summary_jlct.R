summary_jlct <- function(jlct){
  summary <- list(); cg_list <- NULL
  for(l in 1:length(jlct)){
    temp <- jlct[[l]]
    
    p_cg_list <- cg_list
    cg_list <- unlist(sapply(sapply(temp, '[', 'params'), '[', 'gamma'))
    names(cg_list) <- names(temp)
    
    summary_c <- list()
    for(c in 1:length(temp)){
      if(l==1){
        cum_gam <- cg_list[c]
      }else{
        c_node <- names(temp)[c]
        p_node <- paste(unlist(strsplit(c_node, ''))[1:(l+2)], collapse = '')
        cum_gam <- cg_list[c]*p_cg_list[which(names(p_cg_list)==p_node)]
        cg_list[c] <- cum_gam
      }
      
      summary_c[[c]] <- list(var = temp[[c]]$var, gamma = round(temp[[c]]$params$gamma, 3),
                             cum_gamma = round(as.numeric(cum_gam), 3),
                             eta = lapply(temp[[c]]$params$eta, round, 3))
    }
    names(summary_c) <- names(temp)
    summary[[l]] <- summary_c
  }
  names(summary) <- names(jlct)
  
  return(summary)
}