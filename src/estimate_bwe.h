#include <RcppArmadillo.h>


Rcpp::List estimate_bwe(
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
) ;