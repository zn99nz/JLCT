#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::plugins('cpp11')]]
//[[Rcpp::export]]

List post_UD(int njoint, int nitem,
             int n, int nc, NumericVector nclass,
             NumericVector gamma, List eta, List rho,
             List data, NumericVector weight){
  
  List post_jc(n); // √÷¡æ posterior
  NumericVector log_like(n); // likelihood for each subject
  
  NumericMatrix up(njoint, n); // upward prob
  List down(n); // downward prob
  
  
  // UD algorithm
  for(int i=0; i<n; i++){ // for each subject
    
    // (1) calculating up_prob, down_prob
    
    NumericVector up_prod(njoint); // product part of up_prob (1~nc)
    up_prod.fill(1);
    
    List down_temp(nc); // object for ith subject's down_prob
    for(int c=0; c<nc; c++){
      NumericMatrix mat_1(nclass[c], njoint);
      mat_1.fill(1);
      down_temp[c] = mat_1;
    }
    
    for(int c=0; c<nc; c++){
      NumericMatrix eta_c= eta[c];
      NumericMatrix rho_c = rho[c];
      NumericMatrix data_c = data[c];
      
      NumericVector up_sum(njoint);
      up_sum.fill(0);
      
      // object for ith subject's down_prob with cth class variable
      NumericMatrix down_temp_c(nclass[c], njoint); 
      
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
        
        down_temp_c(cl,_) = up_sum_temp; // eta*rho for specific class of cth class variable
      } // for cl-th class for cth class variable
      
      up_prod = up_prod*up_sum; // prod (1~nc) for up_prob
      
      // calculating down_prob
      for(int c1=0; c1<nc; c1++){
        NumericMatrix down_temp_mat = down_temp[c1];
        if(c1 != c){
          for(int cl1=0; cl1<nclass[c1]; cl1++){
            down_temp_mat(cl1,_) = down_temp_mat(cl1,_)*up_sum;
          }
          down_temp[c1] = down_temp_mat;
        }else{
          for(int cl1=0; cl1<nclass[c1]; cl1++){
            down_temp_mat(cl1,_) = down_temp_mat(cl1,_)*down_temp_c(cl1,_);
          }
          down_temp[c1] = down_temp_mat;
        }
      }
      
    } // for cth class variable
    
    // for posterior (up)
    up(_,i) = (up_prod*gamma)/sum(up_prod*gamma); 
    NumericVector up_temp = up(_,i);
    log_like[i] = sum(up_prod*gamma);
    
    // for posterior (down)
    for(int c=0; c<nc; c++){
      NumericMatrix down_temp_post = down_temp[c];
      for(int j=0; j<njoint; j++){
        double down_sum = sum(down_temp_post(_,j));
        if(down_sum == 0){down_sum = 1e-100;}
        down_temp_post(_,j) = down_temp_post(_,j)/down_sum; 
      }
      down_temp[c] = down_temp_post;
    }
    down[i] = down_temp;
    
    
    // (2) calculating ith subject's marginal posterior by up_temp*down_temp
    List post_jc_c(nc);
    
    for(int c=0; c<nc; c++){
      NumericMatrix post_temp = down_temp[c];
      for(int cl=0; cl<nclass[c]; cl++){
        post_temp(cl,_) = up_temp*post_temp(cl,_);
      }
      post_jc_c[c] = post_temp;
    }
    post_jc[i] = post_jc_c;
  } // for ith subject
  
  // result : post_jc, up, weighted log-likelihood
  double w_log_like = sum(log(log_like)*weight);
  
  List result(3);
  result[0] = post_jc; result[1] = up; result[2] = w_log_like;
  
  return result;
  
} // function

