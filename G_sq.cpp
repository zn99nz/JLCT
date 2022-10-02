#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::plugins('cpp11')]]
//[[Rcpp::export]]

NumericVector G_sq_ex(int njoint, int nitem,
             int n, int nc, NumericVector nclass,
             NumericVector gamma, List eta, List rho,
             List data, NumericVector weight){
  
  NumericVector like(n); // likelihood for each subject
  
  NumericMatrix up(njoint, n); // upward prob
  
  // UD algorithm
  for(int i=0; i<n; i++){ // for each subject
    
    // (1) calculating up_prob, down_prob
    
    NumericVector up_prod(njoint); // product part of up_prob (1~nc)
    up_prod.fill(1);
    
    for(int c=0; c<nc; c++){
      NumericMatrix eta_c= eta[c];
      NumericMatrix rho_c = rho[c];
      NumericMatrix data_c = data[c];
      
      NumericVector up_sum(njoint);
      up_sum.fill(0);
      
      for(int cl=0; cl<nclass[c]; cl++){
        // rho part
        double rho_part = 1;
        LogicalVector mis = is_na(data_c(i,_));
        for(int item=0; item<nitem; item++){
          if(mis[item]==FALSE){
            int rsp = data_c(i,item);
            double rho_temp = rho_c(cl,item);
            if(rsp==2){
              rho_temp = 1-rho_temp;
            }
            rho_part = rho_part*rho_temp;
          }
        }
        
        NumericVector up_sum_temp = eta_c(_,cl)*rho_part;
        up_sum = up_sum + up_sum_temp; // summation part of up_prob (for cth class variable)
        
      } // for cl-th class for cth class variable
      
      up_prod = up_prod*up_sum; // prod (1~nc) for up_prob
    } // for cth class variable
    
    // for posterior (up)
    up(_,i) = (up_prod*gamma)/sum(up_prod*gamma); 
    NumericVector up_temp = up(_,i);
    like[i] = sum(up_prod*gamma);
    
  } // for ith subject
  
  return like;
  
} // function

