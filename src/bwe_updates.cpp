// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "etn2.h"
#include "enttn2.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;

// Update z and w, the two truncated normals


arma::mat update_z(
  const arma::mat &beta,
  const arma::mat &x,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ
) {
  arma::mat z(nI, nJ, arma::fill::zeros) ;
//#pragma omp parallel for
  for (int i = 0; i < nI; i++) {
    #pragma omp parallel for
    for (int j = 0; j < nJ; j++) {
      double z_star = 0.0 ;
      arma::vec xi = x.col(i) ;
      arma::vec betaj = beta.col(j) ;
      z_star = as_scalar(xi.t() * betaj);
      double low = mat_Y(i, j) == 1 ? 0.0 : R_NegInf ;
      double high = mat_Y(i, j) == 0 ? 0.0 : R_PosInf ;
      
      z(i, j) = etn2(z_star, 1.0, low, high) ;
      
    }
  }
  return(z) ;
}


arma::mat update_Exx( //note, rename to sum Exx
  arma::mat &var_x,
  arma::mat &x,
  const int nK,
  const int nI){
    arma::mat tmp(nK, nK, arma::fill::zeros) ;
    tmp.fill(nI) ;
    arma::mat Exx = tmp % var_x + x * x.t() ;
    return(Exx) ;
}


arma::mat update_Ebb( // x and beta are Ex and Ebeta
  arma::mat &var_beta,
  arma::mat &beta,
  const int nJ, 
  const int nK){
    arma::mat tmp(nK, nK, arma::fill::zeros) ;
    tmp.fill(nJ) ;
    arma::mat Ebb = tmp % var_beta + beta * beta.t();
    return(Ebb) ;
}


arma::mat update_var_x(
  const arma::mat &xsigma,
  const arma::mat&Ebb
){
  // Calculate A, the variance of X
  arma::mat inv_x_sigma = inv_sympd(xsigma); 
  //#pragma omp parallel for
  arma::mat var_x = inv_sympd( inv_x_sigma + Ebb);
  return(var_x);
}



arma::mat update_x(
  const arma::mat &beta,
  const arma::mat &z,
  const arma::mat &var_x,
  const int nK,
  const int nI,
  const int nJ
){
  arma::mat x(nK, nI, arma::fill::zeros);
  arma::mat tmp(nK, nI, arma::fill::zeros);
  arma::mat inv_var_x = (var_x) ;
  for (int i = 0; i < nI; i++) {
 #pragma omp parallel for
    for (int j = 0; j < nJ; j++) {
        tmp.col(i) += (beta.col(j) * z(i, j));
    }
  }
  
  #pragma omp parallel for
  for(int i = 0; i < nI; i++){
    x.col(i) = inv_var_x * tmp.col(i);
  }
  
  return(x);

}



arma::mat update_var_beta(
  const arma::mat &Exx,
  const arma::mat &betasigma
){
  arma::mat inv_beta_sigma = inv_sympd(betasigma) ; 
  arma::mat var_beta = inv_sympd(inv_beta_sigma + Exx );
  return(var_beta) ;
}



arma::mat update_beta(
  const arma::mat &var_beta,
  const arma::mat &x,
  const arma::mat &z,
  const int nI,
  const int nJ,
  const int nK
){
  arma::mat beta(nK, nJ, arma::fill::zeros);
  arma::mat tmp(nK, nJ, arma::fill::zeros);
  for(int j = 0; j < nJ; j++){
#pragma omp parallel for
    for(int i = 0; i < nI; i++){
    tmp.col(j) += x.col(i) * z(i, j) ;
    }
  }
//  arma::mat inv_var_beta = (var_beta) ;
#pragma omp parallel for
  for(int j = 0; j < nJ; j++){
    beta.col(j) = var_beta * tmp.col(j) ;
  }
  return(beta) ;
}




arma::mat update_z_ns(
    const arma::mat &z,
    const arma::mat &beta,
    const arma::mat &x,
    const arma::sp_mat &mat_Y,
    const int nI,
    const int nJ,
    const arma::mat &p, //unigram dist
    const int prop_ns
) {
  arma::mat z_new = z ;
  for (int i = 0; i < nI; i++) {
    // negative sampling
    arma::mat thisy(mat_Y.row(i));
    int pos_count = sum(thisy.elem( find( thisy == 1 ) )) ;
    int neg_count = prop_ns * pos_count ;
    int samp_count;
    if (neg_count >= nJ - pos_count) {
      samp_count = nJ ;
    } else {
      samp_count = pos_count + neg_count ;
    }
    arma::uvec samp_index(samp_count) ;
    if (neg_count >= nJ - pos_count) {
      for (int n = 0; n < nJ; n++) {
        samp_index[n] = n ;
      }
    } else {
      arma::uvec pos_elements = find( thisy == 1 ) ;
      arma::uvec neg_elements = find( thisy == 0 ) ;
      arma::vec neg_elements_probs = p(neg_elements) ;
      arma::uvec neg_samples = Rcpp::RcppArmadillo::sample_main(
        neg_elements, neg_count, false, neg_elements_probs) ;
      for (int n = 0; n < pos_count; n++) {
        samp_index[n] = pos_elements(n) ;
      }
      for (int n = 0; n < neg_count; n++) {
        samp_index[n + pos_count] = neg_samples(n) ;
      }
    }
#pragma omp parallel for
    for (int j = 0; j < samp_count; j++) {
      double z_star = 0.0 ;
      arma::vec xi = x.col(i) ;
      arma::vec betaj = beta.col(samp_index[j]) ;
      z_star = as_scalar(xi.t() * betaj);
      double low = mat_Y(i, samp_index[j]) == 1 ? 0.0 : R_NegInf ;
      double high = mat_Y(i,samp_index[j]) == 0 ? 0.0 : R_PosInf ;
      z_new(i, samp_index[j]) = etn2(z_star, 1.0, low, high) ;
    }
  }
  return(z_new) ;
}


arma::mat center_x(const arma::mat &x, 
  const int nI, 
  const int nK
  ){
  arma::mat tmp(nK, nI) ;
  double col_avg ;
    for (int i = 0; i < nI; i++) {
      col_avg = mean(x.col(i)) ;
#pragma omp parallel for
        for (int k = 0; k < nK; k++) {
          tmp(k, i) = col_avg ;
        }
    }
    return(tmp) ;
}

// ARD prior for beta

double update_cb(
  const int nJ, 
  const double cb_0
){
  double cb = cb_0 + (nJ * .5) ; 
  return(cb) ;
}


arma::vec update_db(
  const arma::mat &beta,
  const int nK,
  const double db_0
){
  arma::vec d_vals(nK, arma::fill::zeros) ;
  double beta_norm = 0.0 ;
#pragma omp parallel for
  for (int k = 0; k < nK; k++) {
    beta_norm =  norm(beta.row(k), 2) ;
    d_vals(k) =  (db_0 + (pow(beta_norm, 2) * .5)) ;
  }
  return(d_vals) ;
}


arma::mat update_beta_sigma(
  const int nK, 
  const double cb,
  arma::vec db
){
  arma::mat tmp(nK, nK, arma::fill::zeros) ;
#pragma omp parallel for
  for (int k = 0; k < nK; k++) {
    tmp(k, k) =  cb/db(k) ;
  }
  return(tmp) ;
} 
  
// ARD prior for R
double update_cx(
  const int nI, 
  const double cx_0
){
  double cx = cx_0 + (nI * .5) ; 
  return(cx) ;
}

arma::vec update_dx(
  const arma::mat &x,
  const int nK,
  const double dx_0
){
  arma::vec d_vals(nK, arma::fill::zeros) ;
  double x_norm = 0.0 ;
#pragma omp parallel for
  for (int k = 0; k < nK; k++) {
    x_norm =  norm(x.row(k), 2) ;
    d_vals(k) =  (dx_0 + (pow(x_norm, 2) * .5)) ;
  }
  return(d_vals) ;
}


arma::mat update_x_sigma(
  const int nK, 
  const double cx,
  arma::vec dx
){
  arma::mat tmp(nK, nK, arma::fill::zeros) ;
#pragma omp parallel for
  for (int k = 0; k < nK; k++) {
    tmp(k, k) =  cx/dx(k) ;
  }
  return(tmp) ;
} 
  