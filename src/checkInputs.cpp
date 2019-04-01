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
                  unsigned int nJ,
                  unsigned int nI,
                  unsigned int nK,
                  unsigned int threads
                  ) {

    if (maxit < 2) {
        throw std::runtime_error("Max. iterations not > 1.") ;
    }
    if (checkfreq < 1) {
        throw std::runtime_error("Check frequency not positve.") ;
    }

    // if (threads < 0) {
    //     throw std::runtime_error("Number of threads not non-negative.") ;
    // }
    if (nK <= 0) {
        throw std::runtime_error("Number of dimensions not positive.") ;
    }

    // Check Dimensions
    if (verbose) {
        Rcpp::Rcout << "-" << nK << " Dimensional Inputs" << std::endl ;
    }

    //// Priors X
    if ((xmu.n_rows != nI) |
        (xmu.n_cols != nK)
        ) {
        throw std::runtime_error("X prior mean not D x 1.") ;
    }
    if ((xsigma.n_rows != nK) |
        (xsigma.n_cols != nK)
        ) {
        throw std::runtime_error("X prior covariance not D x D.") ;
    }
    //// Priors Alpha, Beta
    if ((betamu.n_rows != nK) |
        (betamu.n_cols != 1)
        ) {
        throw std::runtime_error("Beta prior mean not (D + 1) x 1.") ;
    }
    if ((betasigma.n_rows != nK) |
        (betasigma.n_cols != nK)
        ) {
        throw std::runtime_error("Beta prior covariance not (D + 1) x (D  + 1)") ;
    }

    //// Starts X
    if ((x.n_rows != nI) |
        (x.n_cols != nK)
        ) {
        throw std::runtime_error("X starts not N x D.") ;
    }
    //// Starts Alpha, Beta
    if ((beta.n_rows != nJ) |
        (beta.n_cols != nK)
        ) {
        throw std::runtime_error("Beta starts not J X D.") ;
    }

    // Check Positive-Definiteness of Prior Variances
    arma::mat R ;
    bool test = arma::chol(R, xsigma) ;
    if (!test) {
        throw std::runtime_error("X prior covariance not positive-definite.") ;
    }
    test = arma::chol(R, betasigma) ;
    if (!test) {
        throw std::runtime_error("Beta prior covariance not positive-definite.") ;
    }
}
