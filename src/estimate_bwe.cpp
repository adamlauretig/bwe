// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-
// BWE with ARD prior on beta
#define DEBUG false

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "checkInputs.h"
#include "bwe_updates.h"
#include "elbo_values.h"

using namespace Rcpp ;

// [[Rcpp::export()]]
List estimate_bwe(
  arma::sp_mat mat_Y,
  unsigned int nJ,
  unsigned int nI,
  unsigned int nK,
  arma::mat xmu,
  arma::mat xsigma,
  arma::mat x_init,
  arma::mat betamu,
  arma::mat betasigma,
  arma::mat beta_init,
  bool verbose,
  unsigned int maxit,
  double thresh,
  unsigned int checkfreq,
  unsigned int threads,
  arma::mat p,
  unsigned int prop_ns,
  double cx_0,
  double dx_0,
  double cb_0,
  double db_0
){
  
  
   // Init Qtys
  //// Data Parameters
//  unsigned int nYwords = sum(mat_Y, 1) ;
//  unsigned int nYN = (nD % nJ) - nYwords;

  //Init containers
  
  arma::mat z(nI, nJ, arma::fill::zeros) ;
  arma::vec dx(nK, arma::fill::zeros) ;
  arma::mat Exx(nK, nK, arma::fill::eye) ;
  arma::mat Ebb(nK, nK, arma::fill::eye) ;
  arma::mat var_x(nK, nK, arma::fill::eye) ;
  arma::vec db(nK, arma::fill::zeros) ;
  arma::mat var_beta(nK, nK, arma::fill::eye) ;
  arma::mat x = x_init;
  arma::mat beta = beta_init ;
  
  // "old" containers
  arma::mat oldx = x ;
  arma::mat oldbeta = beta ;
  double old_elbo = 0.0 ;
  // double old_elbo_ns = 0.0 ;
  arma::mat convtrace(maxit, 2, arma::fill::zeros) ;
  arma::vec elbo_tracker(maxit, arma::fill::zeros) ;
  double elbo_est = 0.0 ;
  double Eystar = 0.0 ; 
  
  // c parameter for alpha_b and alpha_x, since it doesn't update
  double cb = 0.0 ;
  cb = update_cb(nJ, cb_0) ;
  double cx = 0.0 ;
  cx = update_cx(nI, cx_0) ;
  
  // arma::vec elbo_tracker_ns(maxit, arma::fill::zeros) ;
  //// Admin
  // unsigned int threadsused = 0 ;
  //bool withstats = false ;
  
  
  unsigned int counter = 0 ;
  bool isconv = false ;
  
#ifdef _OPENMP
  if (threads > 0) {
    omp_set_num_threads( threads );
    // threadsused = omp_get_max_threads() ;
  }
  #endif
  
  while (counter < maxit) {
    counter++ ;
    z = update_z_ns(z, beta, x, mat_Y, nI, nJ, p, prop_ns) ;
    
    var_x = update_var_x(xsigma, Ebb) ;
    dx = update_db(x, nK, dx_0) ;
    xsigma = update_x_sigma(nK, cx, dx) ;
    x = update_x(beta, z, var_x, nK, nI, nJ) ;
    Exx = update_Exx(var_x, x, nK, nI) ;
    db = update_db(beta, nK, db_0) ;
    betasigma = update_beta_sigma(nK, cb, db) ;
    var_beta = update_var_beta(Exx, betasigma) ;
    beta = update_beta(var_beta, x, z, nI, nJ, nK) ;
    Ebb = update_Ebb(var_beta, beta, nJ, nK) ;
    
    
    // convergence
  elbo_est = Eystar + p_z(z, beta, x, Exx, Ebb, mat_Y, nI, nJ) + 
    p_x(Exx, xsigma, nI, nK) + 
    p_beta(Ebb, betasigma, nJ, nK) + 
    p_alpha_x(cx, dx, nK, cx_0, dx_0) +
    p_alpha_b(cb, db, nK, cb_0, db_0) - 
    E_qz(z, beta, x, mat_Y, nI, nJ) - 
    E_qx(var_x, nK, nI) - 
    E_qbeta(var_beta, nK, nJ) -
    E_qalpha_b(cx, dx, nK) - 
    E_qalpha_b(cb, db, nK)
    ;
    
    
    
    
    
// 
//     double elbo_est_ns = elbo_ns( 
//       z, beta, x, mat_Y, var_x, var_beta, Exx, Ebb, xsigma, betasigma, 
//       nI, nJ, nK, p, prop_ns) ;

       
    double delta_elbo =  (elbo_est - old_elbo)/std::abs(old_elbo) ;
    // double delta_elbo_ns = elbo_est_ns - old_elbo_ns;
    // convtrace(counter - 1, 0) = devx ;
    // convtrace(counter - 1, 1) = devbeta ;
    convtrace(counter - 1, 0) = rank(x) ;
    convtrace(counter - 1, 1) = delta_elbo;
    // convtrace(counter - 1, 2) = delta_elbo_ns;

    elbo_tracker(counter - 1) = elbo_est ;  
    // elbo_tracker_ns(counter - 1) = elbo_est_ns ; 
    isconv = ((delta_elbo) < thresh) ;

    // if(prop_ns > 0){
    //   isconv = (std::abs(delta_elbo_ns) < thresh) ;
    // } else{
    // }
    if (counter % checkfreq == 0 & counter > 1) {
      R_CheckUserInterrupt() ;
      if (verbose) {
        Rcout << "Iteration: " << counter ;
        Rcout << "     conv: " << convtrace.row(counter - 1) ;
      }
    }

    if ((counter > 2) & isconv) {
      break ;
    } else {
      oldx = x ;
      oldbeta = beta ;
      old_elbo = elbo_est  ;
      // old_elbo_ns = elbo_est_ns ;
    }
  }

  List ret ;
  // ret["mat_Y"] = mat_Y ;
  ret["x"] = x ;
  ret["beta"] = beta ;  
  // ret["nJ"] = nJ ;
  // ret["nI"] = nI ;
  // ret["z"] = z ;
  ret["var_x"] = var_x ;
  ret["xsigma"] = xsigma ;
  ret["cx"] = cx ;
  ret["dx"] = dx ;
  ret["cb"] = cb ;
  ret["db"] = db ;
  ret["var_beta"] = var_beta ;
  ret["betasigma"] = betasigma ;
  // ret["Exx"] = Exx ;
  // ret["Ebb"] = Ebb ;
  ret["convergence"] = convtrace ;
  ret["elbo"] = elbo_tracker ;
  // ret["elbo_ns"] = elbo_tracker_ns ;
  return(ret); 
}