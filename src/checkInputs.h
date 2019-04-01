// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

void checkInputs (const arma::sp_mat &mat_Y,
                  const arma::mat &beta,
                  const arma::mat &x,
                  const arma::mat &xmu,
                  const arma::mat &xsigma,
                  const arma::mat &betamu,
                  const arma::mat &betasigma,
                  bool verbose,
                  unsigned int maxit,
                  double thresh,
                  unsigned int checkfreq,
                  unsigned int threads,
                  unsigned int nJ,
                  unsigned int nI,
                  unsigned int nK
                  ) ;
