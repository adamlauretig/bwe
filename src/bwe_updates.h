// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-


#include <RcppArmadillo.h>

arma::mat update_z(
  const arma::mat &beta,
  const arma::mat &x,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ
) ;
  
arma::mat update_Exx(
  arma::mat &var_x,
  arma::mat &x,
  const int nK,
  const int nI
) ;

arma::mat update_Ebb(
  arma::mat &var_beta,
  arma::mat &beta,
  const int nJ, 
  const int nK
) ;

arma::mat update_var_x(
  const arma::mat &xsigma,
  const arma::mat&Ebb
) ;


arma::mat update_x(
  const arma::mat &beta,
  const arma::mat &z,
  const arma::mat &var_x,
  const int nK,
  const int nI,
  const int nJ
) ;

arma::mat update_var_beta(
  const arma::mat &Exx,
  const arma::mat &betasigma
) ;

arma::mat update_beta(
  const arma::mat &var_beta,
  const arma::mat &x,
  const arma::mat &z,
  const int nI,
  const int nJ,
  const int nK
) ;

arma::mat update_z_ns(
  const arma::mat &z,
  const arma::mat &beta,
  const arma::mat &x,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ,
  const arma::mat &p, //unigram dist
  const int prop_ns
) ;

arma::mat update_z_ns_tri(
  const arma::mat &z,
  const arma::mat &beta,
  const arma::mat &x,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ,
  const arma::mat &p, //unigram dist
  const int prop_ns
) ;

double elbo(
    const arma::mat &z,
    const arma::mat &beta,
    const arma::mat &x,
    const arma::sp_mat &mat_Y,
    const arma::mat &var_x,
    const arma::mat &var_beta,
    const arma::mat &Exx,
    const arma::mat &Ebb,
    const arma::mat &xsigma,
    const arma::mat &betasigma,
    const int nI,
    const int nJ,
    const int nK
) ;

double elbo_ns(
    const arma::mat &z,
    const arma::mat &beta,
    const arma::mat &x,
    const arma::sp_mat &mat_Y,
    const arma::mat &var_x,
    const arma::mat &var_beta,
    const arma::mat &Exx,
    const arma::mat &Ebb,
    const arma::mat &xsigma,
    const arma::mat &betasigma,
    const int nI,
    const int nJ,
    const int nK,    
    const arma::mat &p, //unigram dist
    const int prop_ns
) ;

arma::mat center_x(const arma::mat &x, 
  const int nI, 
  const int nK) ;

double update_cb(
  const int nJ, 
  const double cb_0
) ;

arma::vec update_db(
  const arma::mat &beta,
  const int nK,
  const double db_0
) ;

arma::mat update_beta_sigma(
  const int nK, 
  const double cb,
  arma::vec db
) ;

double update_cx(
  const int nI, 
  const double cx_0
) ;

arma::vec update_dx(
  const arma::mat &x,
  const int nK,
  const double dx_0
) ;

arma::mat update_x_sigma(
  const int nK, 
  const double cx,
  arma::vec dx
) ;
