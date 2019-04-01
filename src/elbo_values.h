// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-


#include <RcppArmadillo.h>

double E_qz(
    const arma::mat &z,
    const arma::mat &beta,
    const arma::mat &x,
    const arma::sp_mat &mat_Y,
    const int nI,
    const int nJ
) ;

double E_qx(
  const arma::mat &var_x,
  const int nK,
  const int nI
) ;

double E_qbeta(
  const arma::mat &var_beta,
  const int nK,
  const int nJ
) ;

double E_qalpha_x(
  const double cx,
  arma::vec dx,
  const int nK
) ;

double E_qalpha_b(
  const double cb,
  arma::vec db,
  const int nK
) ;

double p_z(
  const arma::mat &z,
  const arma::mat &beta,
  const arma::mat &x,
  const arma::mat &Exx,
  const arma::mat &Ebb,
  const arma::sp_mat &mat_Y,
  const int nI,
  const int nJ
) ;

double p_x(
  const arma::mat &Exx,
  const arma::mat &xsigma,
  const int nI,
  const int nK
) ;

double p_beta(
  const arma::mat &Ebb,
  const arma::mat &betasigma,
  const int nJ,
  const int nK
) ;

double p_alpha_b(
  const double cb,
  arma::vec db,
  const int nK,
  const double cb_0,
  const double db_0
) ;

double p_alpha_x(
  const double cx,
  arma::vec dx,
  const int nK,
  const double cx_0,
  const double dx_0
) ;