// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "etn2.h"
#include "enttn2.h"

using namespace Rcpp ;
using namespace RcppArmadillo ;


// qz entropy

double E_qz(
    const arma::mat &z,
    const arma::mat &beta,
    const arma::mat &x,
    const arma::sp_mat &mat_Y,
    const int nI,
    const int nJ
){
  double qz = 0.0;
   for (int i = 0; i < nI; i++) {
#pragma omp parallel for reduction(+:qz)
    for (int j = 0; j < nJ; j++) {
      double z_star = 0.0 ;
      arma::vec xi = x.col(i) ;
      arma::vec betaj = beta.col(j) ;
      z_star = as_scalar(xi.t() * betaj);
      double low = mat_Y(i, j) == 1 ? 0.0 : R_NegInf ;
      double high = mat_Y(i, j) == 0 ? 0.0 : R_PosInf ;
      double tmp = enttn2(z_star, 1.0, low, high) ;
      qz += tmp ;
    }
  }
   return(qz) ;
}


// qx entropy

double E_qx(
  const arma::mat &var_x,
  const int nK,
  const int nI
){
  double qx = 0.0 ;
  double tmp = nK * .5 ;
  double tmp2 = std::log(2 * arma::datum::pi * exp(1)) ;
  double tmp3 = .5 * std::log(det(var_x)) ;
  // qx = -nI * (  + (.5 * std::log(det(var_x))));
  
  qx = (tmp * tmp2) + tmp3 ;
  qx = nI * qx ;
  return(qx) ;
}


double E_qbeta(
  const arma::mat &var_beta,
  const int nK,
  const int nJ
){
  double qbeta = 0.0 ;
  double tmp = nK * .5 ;
  double tmp2 = std::log(2 * arma::datum::pi * exp(1)) ;
  double tmp3 = .5 * std::log(det(var_beta)) ;

  qbeta = (tmp * tmp2) + tmp3 ;
  qbeta = nJ * qbeta ;
  return(qbeta) ;
}

//E(q_alpha)

double E_qalpha_b(
  const double cb,
  arma::vec db,
  const int nK
){
  double q_alpha = 0.0 ;
  double q1 = R::lgammafn(cb) ;
  double q2 = R::digamma(cb) ;
  
// #pragma omp parallel for reduction(+:q_alpha)
    for (int k = 0; k < nK; k++) {  
      q_alpha += q1 - ((cb - 1) * q2) - log(db(k)) + cb ;
    }
    return(q_alpha) ;
}


double E_qalpha_x(
  const double cx,
  arma::vec dx,
  const int nK
){
  double q_alpha = 0.0 ;
  double q1 = R::lgammafn(cx) ;
  double q2 = R::digamma(cx) ;
  
// #pragma omp parallel for reduction(+:q_alpha)
    for (int k = 0; k < nK; k++) {  
      q_alpha += q1 - ((cx - 1) * q2) - log(dx(k)) + cx ;
    }
    return(q_alpha) ;
}

//E(p_z)

double p_z(
  const arma::mat &z,
  const arma::mat &beta,
  const arma::mat &x,
  const arma::mat &Exx,
  const arma::mat &Ebb,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ
){
  double pz = 0.0 ;
  for (int i = 0; i < nI; i++) {
  arma::vec E_xi = x.col(i) ;
    
#pragma omp parallel for reduction(+:pz)
    for (int j = 0; j < nJ; j++) {
      arma::vec E_betaj = beta.col(j) ;
      
      double q1 = -(.5) * log( 2.0 * arma::datum::pi) ;
      double q2 = as_scalar(
        (.5) * (pow(z(i, j), 2) - 2 * z(i, j) * E_xi.t() * E_betaj ));
      pz += q1 - q2 ;
    }
  }
  pz = pz + trace(Exx * Ebb) ;
  return(pz) ;
}

// E(p_x)

double p_x(
  const arma::mat &Exx,
  const arma::mat &xsigma,
  const int nI,
  const int nK
){
  double tmp1 = -1 * (nK * .5 * log(2 * arma::datum::pi)) ;
  double tmp2 = .5 * std::log(det(xsigma)) ;
  double tmp3 = nI * (tmp1 - tmp2) ;
  double px = tmp3 - (.5 * trace(xsigma * Exx.t())) ;
  return(px) ;
}

double p_beta(
  const arma::mat &Ebb,
  const arma::mat &betasigma,
  const int nJ,
  const int nK
){
  double tmp1 = -1 * (nK * .5 * log(2 * arma::datum::pi)) ;
  double tmp2 = .5 * std::log(det(betasigma)) ;
  double tmp3 = nJ * (tmp1 - tmp2) ;
  double pbeta = tmp3 - (.5 * trace(betasigma * Ebb.t())) ;
  return(pbeta) ;
}


double p_alpha_b(
  const double cb,
  arma::vec db,
  const int nK,
  const double cb_0,
  const double db_0
){
  double p_alpha_val = 0.0 ;
  double q1 = (cb_0 * log(db_0)) ;
  double q2 = R::lgammafn(cb) ;
// #pragma omp parallel for reduction(+:q_alpha)
    for (int k = 0; k < nK; k++) {  
      p_alpha_val += q1 + (cb_0 - 1) * (
        R::digamma(cb) - log(db(k))) - (db_0 * (cb/db(k))) + q2 ;
    }
  return (p_alpha_val) ;
}

double p_alpha_x(
  const double cx,
  arma::vec dx,
  const int nK,
  const double cx_0,
  const double dx_0
){
  double p_alpha_val = 0.0 ;
  double q1 = (cx_0 * log(dx_0)) ;
  double q2 = R::lgammafn(cx) ;
// #pragma omp parallel for reduction(+:q_alpha)
    for (int k = 0; k < nK; k++) {  
      p_alpha_val += q1 + (cx_0 - 1) * (
        R::digamma(cx) - log(dx(k))) - (dx_0 * (cx/dx(k))) + q2 ;
    }
  return (p_alpha_val) ;
}

