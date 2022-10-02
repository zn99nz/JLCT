rm(list = ls())

library(haven)
yrbs2015 <- read_sas('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\YRBS_2015\\yrbs2015.sas7bdat')
yrbs2015 <- as.data.frame(yrbs2015)
yrbs2015

# 1 : age 
# 2 : sex

# violent behavior : 13, 16, 17, 18, 19, 
# depression : 26~30
# drug : 31, 33, 41, 43, 47
# sexual behavior : 60~64

# 필요한 변수만 선택
rq_col <- c(1:2, 13, 16:19, 26:30, 31, 33, 41, 43, 47, 60:64)
yrbs_1 <- yrbs2015[,rq_col]
for(i in 1:ncol(yrbs_1)){
  yrbs_1[which(yrbs_1[,i] == ''),i] <- NA
  yrbs_1[,i] <- as.numeric(yrbs_1[,i])
}

# 나이, 성별 모두 답변한 학생만
yrbs_2 <- yrbs_1[complete.cases(yrbs_1[,1:2]),]
dim(yrbs_2)

# 16~18세 남학생만
cond <- (yrbs_2$Q1 %in% 5:7 & yrbs_2$Q2 == 2)
yrbs_3 <- yrbs_2[cond,]
dim(yrbs_3) # n = 5058

yrbs_4 <- yrbs_3[,1:2]

resp_trans <- function(x, resp1, resp2){
  # 1 : yes
  # 2 : no

  resp_new <- rep(NA, length(x))
  resp_new[which(x %in% resp1&!is.na(x))] <- 1
  resp_new[which(x %in% resp2&!is.na(x))] <- 2
  
  return(resp_new)
}


#### violent behavior
viol <- which(colnames(yrbs_3) %in% paste0('Q', c(13,16:19)))
yrbs_viol <- yrbs_3[,viol]

yrbs_4[,viol[1]] <- resp_trans(yrbs_viol[,1], 2:5, 1)
yrbs_4[,viol[2]] <- resp_trans(yrbs_viol[,2], 2:5, 1)
yrbs_4[,viol[3]] <- resp_trans(yrbs_viol[,3], 2:8, 1)
yrbs_4[,viol[4]] <- resp_trans(yrbs_viol[,4], 2:8, 1)
yrbs_4[,viol[5]] <- resp_trans(yrbs_viol[,5], 2:5, 1)

colnames(yrbs_4)[viol] <- paste0('viol', 1:5)
viol_NA <- which(apply(is.na(yrbs_4[,viol]), 1, any))


#### depression
dprs <- which(colnames(yrbs_3) %in% paste0('Q', c(26:30)))
yrbs_dprs <- yrbs_3[,dprs]

yrbs_4[,dprs[1]] <- resp_trans(yrbs_dprs[,1], 1, 2)
yrbs_4[,dprs[2]] <- resp_trans(yrbs_dprs[,2], 1, 2)
yrbs_4[,dprs[3]] <- resp_trans(yrbs_dprs[,3], 1, 2)
yrbs_4[,dprs[4]] <- resp_trans(yrbs_dprs[,4], 2:5, 1)
yrbs_4[,dprs[5]] <- resp_trans(yrbs_dprs[,5], 2, c(1,3))

colnames(yrbs_4)[dprs] <- paste0('dprs', 1:5)
dprs_NA <- which(apply(is.na(yrbs_4[,dprs]), 1, any))


#### drug
drug <- which(colnames(yrbs_3) %in% paste0('Q', c(31, 33, 41, 43, 47)))
yrbs_drug <- yrbs_3[,drug]

yrbs_4[,drug[1]] <- resp_trans(yrbs_drug[,1], 1, 2)
yrbs_4[,drug[2]] <- resp_trans(yrbs_drug[,2], 2:7, 1)
yrbs_4[,drug[3]] <- resp_trans(yrbs_drug[,3], 2:7, 1)
yrbs_4[,drug[4]] <- resp_trans(yrbs_drug[,4], 2:7, 1)
yrbs_4[,drug[5]] <- resp_trans(yrbs_drug[,5], 2:7, 1)

colnames(yrbs_4)[drug] <- paste0('drug', 1:5)
drug_NA <- which(apply(is.na(yrbs_4[,drug]), 1, any))


#### sexual behavior
sbhv <- which(colnames(yrbs_3) %in% paste0('Q', c(60:64)))
yrbs_sbhv <- yrbs_3[,sbhv]

yrbs_4[,sbhv[1]] <- resp_trans(yrbs_sbhv[,1], 1, 2)
yrbs_4[,sbhv[2]] <- resp_trans(yrbs_sbhv[,2], 2:3, c(1, 4:8))
yrbs_4[,sbhv[3]] <- resp_trans(yrbs_sbhv[,3], 5:7, 1:4)
yrbs_4[,sbhv[4]] <- resp_trans(yrbs_sbhv[,4], 3:8, 1:2)
yrbs_4[,sbhv[5]] <- resp_trans(yrbs_sbhv[,5], 2, c(1,3))

colnames(yrbs_4)[sbhv] <- paste0('sbhv', 1:5)
sbhv_NA <- which(apply(is.na(yrbs_4[,sbhv]), 1, any))

ind_NA <- unique(c(viol_NA, dprs_NA, drug_NA, sbhv_NA))
yrbs_5 <- yrbs_4[-ind_NA,]
nrow(yrbs_5)




yrbs_data <- list(); length(yrbs_data) <- 4
names(yrbs_data) <- c('violent_behavior', 'depression', 
                      'drug', 'sexual_behavior')

ind_list <- matrix(c(viol, dprs, drug, sbhv), 5,4)

for(i in 1:4){
  yrbs_data[[i]] <- as.matrix(yrbs_5[,ind_list[,i]])
}

save(yrbs_data, file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\yrbs_data.RData')

#####################
####     LCA     ####
#####################
n=nrow(yrbs_data[[1]])

library(poLCA)
source('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\data generate.r')

#### 1. violent behavior
f_viol <- cbind(viol1, viol2, viol3, viol4, viol5) ~ 1
viol_list <- list()

viol_mat <- matrix(0, 4, 3); colnames(viol_mat) <- c('AIC', 'BIC', 'p_val')
maxiter <- c(5000, 5000, 10000, 10000)
for(nc in 2:5){
  viol_lca <- poLCA(f_viol, as.data.frame(yrbs_data[[1]]), 
                    nclass = nc, maxiter = maxiter[nc-1], tol = 1e-6, nrep = 10, verbose = FALSE)
  viol_list[[nc-1]] <- viol_lca
  print(paste('viol_nc :', nc))
  
  # AIC, BIC
  viol_mat[(nc-1),1] <- viol_lca$aic; viol_mat[(nc-1),2] <- viol_lca$bic
  
  # Bootstrap Gsq
  Gsq <- rep(0, 100)
  for(i in 1:100){
    data_boot <- dgen_LCA(viol_lca$P, viol_lca$probs, n, colnames(yrbs_data[[1]]))
    viol_boot <- poLCA(f_viol, data_boot, nclass = nc, tol = 1e-6, verbose = FALSE)
    Gsq[i] <- viol_boot$Gsq
    if(i%%10 == 0){print(paste('   boot :', i))}
  }
  viol_mat[(nc-1),3] <- mean(Gsq > viol_lca$Gsq)
}
viol_mat

#### 2. depression
f_dprs <- cbind(dprs1, dprs2, dprs3, dprs4, dprs5) ~ 1
dprs_list <- list()

dprs_mat <- matrix(0, 4, 3); colnames(dprs_mat) <- c('AIC', 'BIC', 'p_val')
maxiter <- c(5000, 5000, 5000, 5000)
for(nc in 2:5){
  dprs_lca <- poLCA(f_dprs, as.data.frame(yrbs_data[[2]]), 
                    nclass = nc, maxiter = maxiter[nc-1], tol = 1e-6, nrep = 20, verbose = FALSE)
  
  dprs_list[[nc-1]] <- dprs_lca
  print(paste('dprs_nc :', nc))
  
  # AIC, BIC
  dprs_mat[(nc-1),1] <- dprs_lca$aic; dprs_mat[(nc-1),2] <- dprs_lca$bic
  
  # Bootstrap Gsq
  Gsq <- rep(0, 100)
  for(i in 1:100){
    data_boot <- dgen_LCA(dprs_lca$P, dprs_lca$probs, n, colnames(yrbs_data[[2]]))
    dprs_boot <- poLCA(f_dprs, data_boot, nclass = nc, tol = 1e-6, verbose = FALSE)
    Gsq[i] <- dprs_boot$Gsq
    if(i%%10 == 0){print(paste('   boot :', i))}
  }
  dprs_mat[(nc-1),3] <- mean(Gsq > dprs_lca$Gsq)
}
dprs_mat


#### 3. drug
f_drug <- cbind(drug1, drug2, drug3, drug4, drug5) ~ 1
drug_list <- list()

drug_mat <- matrix(0, 4, 3); colnames(drug_mat) <- c('AIC', 'BIC', 'p_val')
maxiter <- c(5000, 5000, 10000, 5000)

for(nc in 2:5){
  drug_lca <- poLCA(f_drug, as.data.frame(yrbs_data[[3]]), 
                    nclass = nc, maxiter = maxiter[nc-1], tol = 1e-6, nrep = 10, verbose = FALSE)
  drug_list[[nc-1]] <- drug_lca
  print(paste('drug_nc :', nc))
  
  # AIC, BIC
  drug_mat[(nc-1),1] <- drug_lca$aic; drug_mat[(nc-1),2] <- drug_lca$bic
  
  # Bootstrap Gsq
  Gsq <- rep(0, 100)
  for(i in 1:100){
    data_boot <- dgen_LCA(drug_lca$P, drug_lca$probs, n, colnames(yrbs_data[[3]]))
    drug_boot <- poLCA(f_drug, data_boot, nclass = nc, tol = 1e-6, verbose = FALSE)
    Gsq[i] <- drug_boot$Gsq
    if(i%%10 == 0){print(paste('   boot :', i))}
  }
  drug_mat[(nc-1),3] <- mean(Gsq > drug_lca$Gsq)
}
drug_mat


#### 4. sexual behavior
f_sbhv <- cbind(sbhv1, sbhv2, sbhv3, sbhv4, sbhv5) ~ 1
sbhv_list <- list()

sbhv_mat <- matrix(0, 4, 3); colnames(sbhv_mat) <- c('AIC', 'BIC', 'p_val')
maxiter <- c(5000, 5000, 10000, 10000)

for(nc in 2:5){nc=4
  sbhv_lca <- poLCA(f_sbhv, as.data.frame(yrbs_data[[4]]), 
                    nclass = nc, maxiter = maxiter[nc-1], tol = 1e-6, nrep = 10, verbose = FALSE)
  sbhv_lca$aic; sbhv_lca$bic; sbhv_lca$Gsq
  
  sbhv_list[[nc-1]] <- sbhv_lca
  print(paste('sbhv_nc :', nc))
  
  # AIC, BIC
  sbhv_mat[(nc-1),1] <- sbhv_lca$aic; sbhv_mat[(nc-1),2] <- sbhv_lca$bic
  
  # Bootstrap Gsq
  Gsq <- rep(0, 100)
  for(i in 1:100){
    data_boot <- dgen_LCA(sbhv_lca$P, sbhv_lca$probs, n, colnames(yrbs_data[[4]]))
    sbhv_boot <- poLCA(f_sbhv, data_boot, nclass = nc, tol = 1e-6, verbose = FALSE)
    Gsq[i] <- sbhv_boot$Gsq
    if(i%%10 == 0){print(paste('   boot :', i))}
  }
  sbhv_mat[(nc-1),3] <- mean(Gsq > sbhv_lca$Gsq)
}
sbhv_mat


rho_list <- list()
# rho : list ((nclass*nitem) * nc)
viol_rho <- matrix(0, 4, 5)
viol_lca <- viol_list[[3]]
for(i in 1:5){viol_rho[,i] <- viol_lca$probs[[i]][,1]}
rho_list[[1]] <- viol_rho[4:1,]

dprs_rho <- matrix(0, 3, 5)
dprs_lca <- dprs_list[[2]]
for(i in 1:5){dprs_rho[,i] <- dprs_lca$probs[[i]][,1]}
rho_list[[2]] <- dprs_rho[3:1,]

drug_rho <- matrix(0, 4, 5)
drug_lca <- drug_list[[3]]
for(i in 1:5){drug_rho[,i] <- drug_lca$probs[[i]][,1]}
rho_list[[3]] <- drug_rho[c(2,3,1,4),]

sbhv_rho <- matrix(0, 3, 5)
sbhv_lca <- sbhv_list[[2]]
for(i in 1:5){sbhv_rho[,i] <- sbhv_lca$probs[[i]][,1]}
rho_list[[4]] <- sbhv_rho[c(3,1,2),]

### LCA result save
fixed_lca <- list()
fixed_lca[[1]] <- viol_list[[3]]
fixed_lca[[2]] <- dprs_list[[2]]
fixed_lca[[3]] <- drug_list[[3]]
fixed_lca[[4]] <- sbhv_list[[2]]
names(fixed_lca) <- names(yrbs_data)
save(fixed_lca, file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\fixed_lca.RData')

lca_diag <- list()
lca_diag[[1]] <- viol_mat
lca_diag[[2]] <- dprs_mat
lca_diag[[3]] <- drug_mat
lca_diag[[4]] <- sbhv_mat
names(lca_diag) <- names(yrbs_data)
save(lca_diag, file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\lca_diag.RData')

### rho result save
save(rho_list, file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\rho_list.RData')
##############################################################################################################
rm(list=ls())

load('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\yrbs_data.RData')
resp_check <- function(x){
  tab <- apply(x, 2, table, useNA = 'ifany')
  tab_ratio <- apply(tab, 2, function(x) return(x/sum(x)))
  return(round(tab_ratio, 3))
}

lapply(yrbs_data, resp_check)

######################
####     JLCT     ####
######################
n=nrow(yrbs_data[[1]])
load('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\fixed_result\\rho_list.RData')
source('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\JLCT_code.r')

jlct_result_newnew <- list()
for(i in 1:5){
  t1 <- Sys.time()
  yrbs_result <- JLCT(yrbs_data, n=nrow(yrbs_data[[1]]), nc=4,
                        nclass=c(4,3,4,3), nitem=5, rho_fix = TRUE, rho_init = rho_list)
  t2 <- Sys.time()
  jlct_result_newnew[[i]] <- yrbs_result
  print(t2-t1)
}

save(jlct_result_newnew, 
     file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\wamp_result\\jlct_result_newnew.RData')
length(jlct_result_newnew)

t1 <- Sys.time()
jlct_result_new2 <- JLCT(yrbs_data, n=nrow(yrbs_data[[1]]), nc=4,
                    nclass=c(4,3,4,3), nitem=5, rho_fix = TRUE, rho_init = rho_list)
t2 <- Sys.time()
save(jlct_result_new2, file = 'C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\wamp_result\\jlct_result_new2.RData')


source('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\summary_jlct.r')

load('C:\\Users\\uccoo\\Desktop\\학교\\study\\JLCT\\wamp_result\\jlct_result_newnew.RData')

jlct_1 <- jlct_result_newnew[[1]]
lapply(jlct_1, names)

summary_jlct(jlct_1)
summary_jlct(jlct_2)

jlct_2 <- jlct_result_newnew[[2]]
summary_jlct(jlct_2)

summary_jlct(jlct_result_new2)
lapply(jlct_result_new2, names)




#########################################################################3

rm(list=ls())




