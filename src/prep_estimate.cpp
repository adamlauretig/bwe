
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "estimate_bwe.h"


// [[Rcpp::export()]]
RcppExport SEXP estimate(SEXP mat_Y_r,
                                 SEXP nJ_r,
                                 SEXP nI_r,
                                 SEXP nK_r,
                                 SEXP xmu_r,
                                 SEXP xsigma_r,
                                 SEXP x_init_r,
                                 SEXP betamu_r,
                                 SEXP betasigma_r,
                                 SEXP beta_init_r,
                                 SEXP verbose_r,
                                 SEXP maxit_r,
                                 SEXP thresh_r,
                                 SEXP checkfreq_r,
                                 SEXP threads_r,
                                 SEXP p_r,
                                 SEXP prop_ns_r,
                                 SEXP cx_0_r,
                                 SEXP dx_0_r,
                                 SEXP cb_0_r,
                                 SEXP db_0_r
                                 ) {
    BEGIN_RCPP
    SEXP resultSEXP ;
  {
    Rcpp::RNGScope __rngScope ;
    Rcpp::traits::input_parameter<arma::sp_mat>::type mat_Y(mat_Y_r) ;
    Rcpp::traits::input_parameter<int>::type nJ(nJ_r) ;
    Rcpp::traits::input_parameter<int>::type nI(nI_r) ;
    Rcpp::traits::input_parameter<int>::type nK(nK_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type xmu(xmu_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type xsigma(xsigma_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type x_init(x_init_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type betamu(betamu_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type betasigma(betasigma_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type beta_init(beta_init_r) ;
    Rcpp::traits::input_parameter<bool>::type verbose(verbose_r) ;
    Rcpp::traits::input_parameter<int>::type maxit(maxit_r) ;
    Rcpp::traits::input_parameter<double>::type thresh(thresh_r) ;
    Rcpp::traits::input_parameter<int>::type checkfreq(checkfreq_r) ;
    Rcpp::traits::input_parameter<int>::type threads(threads_r) ;
    Rcpp::traits::input_parameter<arma::mat>::type p(p_r) ;
    Rcpp::traits::input_parameter<int>::type prop_ns(prop_ns_r) ;
    Rcpp::traits::input_parameter<double>::type cx_0(cx_0_r) ;
    Rcpp::traits::input_parameter<double>::type dx_0(dx_0_r) ;
    Rcpp::traits::input_parameter<double>::type cb_0(cb_0_r) ;
    Rcpp::traits::input_parameter<double>::type db_0(db_0_r) ;

    Rcpp::List result = estimate_bwe(mat_Y,
                                 nJ,
                                 nI,
                                 nK,
                                 xmu,
                                 xsigma,
                                 x_init,
                                 betamu,
                                 betasigma,
                                 beta_init,
                                 verbose,
                                 maxit,
                                 thresh,
                                 checkfreq,
                                 threads,
                                 p,
                                 prop_ns,
                                 cx_0,
                                 dx_0,
                                 cb_0,
                                 db_0
                                 ) ;
    PROTECT(resultSEXP = Rcpp::wrap(result)) ;
  }
  UNPROTECT(1);
  return(resultSEXP) ;
  END_RCPP
    }
