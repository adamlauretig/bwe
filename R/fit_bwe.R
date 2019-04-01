#' Estimate BWE
#' 
#' Takes a collection of documents, converts them to a format for statistical estimation, and then estimates a Bayesian word embedding model using variational inference. 
#' 
#' @param dat A \code{data.table} or \code{data.frame}, where every row is a document.
#' @param text_field Character string. The name of the column in \code{dat} which contains the text.
#' @param context_window Integer. The size of the word co-occurrence window. Defaults to 3.
#' @param minimum_word_count Integer. The minimum word frequency to be included in the corpus. Defaults to 5, but with larger corpora, increase this. 
#' @param character_min Integer. What is the smallest sequence of characters to be considered? Anything smaller than this is discarded. Defaults to 0.
#' @param remove_numbers Logical. Argument to \code{quanteda} to remove numbers from the co-occurrence matrix. Defaults to \code{FALSE}.
#' @param stops Vector. A vector of stopwords to remove from the co-occurrence matrix. Defaults to \code{NULL}, with no words removed.
#' @param K Integer. The number of embedding dimensions. Defaults to five, but should usually be larger.
#' @param start_type Character. How the \code{x} and \code{beta} matrices are initialized. \code{zeroes} initialized with a matrix of zeroes, and \code{random} intializes with draws from a standard normal distribution. \code{svd} calculates the singular value decomposition of the co-occurrence matrix using \code{\link[irlba]{irlba}}, and initializes \code{x} and \code{beta} with the \code{u} and \code{v} values output from \code{irlba}. We recommend using this initializing, as it provides more useful initializations for the model.
#' @param verbose Boolean. Should model convergence be displayed as the model runs? Defaults to \code{TRUE}.
#' @param maxit Integer. Maximum number of iterations to run the model. Default is 5000.
#' @param thresh Numeric. Convergence threshold (Change in ELBO). Defaults to .01.
#' @param checkfreq Integer. How often should verbose output be displayed?
#' @param threads Integer. Number of cores to use, default is 1, threading done via OpenMP.
#' @param prop_ns Integer. Number of negative samples for every positive sample. Default is 5, but this should be adjusted depending on your dataset.
#' @param cx Positive numeric. Hyperparameters for \code{x} covariance. 
#' @param dx Positive numeric. Hyperparameters for \code{x} covariance. 
#' @param cb Positive numeric. Hyperparameters for \code{beta} covariance. 
#' @param db Positive numeric. Hyperparameters for \code{beta} covariance. 
#' @import data.table
#' @import Matrix
#' @import quanteda
#' @import Rcpp
#' @import RcppArmadillo
#' @return A list of objects
#' 
#' \item{data_prep}{A list of the inputs into the \code{bwe} function. Includes sparse co-occurrence matrix, word frequencies, transformed probabilities for negative sampling, and initialization values.}
#' \item{x}{The matrix of word embeddings. Rows are dimensions, columns are words.}
#' \item{beta}{The matrix of context embeddings. Rows are dimensions, columns are words.}
#' \item{var_x}{The covariance matrix for \code{x}.}
#' \item{var_beta}{The covariance matrix for \code{beta}.}
#' \item{xsigma}{The prior, learned from the data, for x.}
#' \item{betasigma}{The prior, learned from the data, for beta.}
#' \item{cx, dx, cb, bd}{The priors on \code{xsigma} and \code{betasigma}, learned from the data.}
#' \item{convergence}{A two column matrix. Column 1 is the rank of \code{x}, column two is the change in the elbo.}
#' \item{elbo}{A vector. The raw elbo values.}
#' 
#' @useDynLib bwe 
#' 
#' 
#' @export

fit_bwe <- function(dat, 
  text_field = NA, 
  context_window = 3L,
  minimum_word_count = 5L,
  K = 5, 
  character_min = 0L,
  remove_numbers = FALSE,
  stops = NULL,
  start_type = c("svd", "zeros", "random"),
  verbose = TRUE,
  maxit = 5000,
  thresh = .01,
  checkfreq = 5,
  prop_ns = 5,
  cx_0 = .01, 
  dx_0 = .01,
  cb_0 = .01,
  db_0 = .01,
  seed = NA){
  
  # context_window <- as.integer(context_window)
  # minimum_word_count <- as.integer(minimum_word_count)
  # K <- as.integer(K)
  # maxit <- as.integer(maxit)
  # checkfreq <- as.integer(checkfreq)
  # 
  # 
  if(is.na(seed)){
    s <- sample.int(n = .Machine$integer.max, size = 1)
  } else{
    s <- seed
  }
  set.seed(s)
  
  d <- data_prep(dat = dat, 
    text_field = text_field, 
    window_to_use = context_window, 
    minimum_word_count = minimum_word_count,
    start_type = start_type, K = K, char_min = character_min, 
    rm_nums = remove_numbers, stops_list = stops)
 
  cat("Initialization complete, starting estimation.\n")
   
  bwe_out <- estimate(
  mat_Y_r = d$mat_Y,
  nJ_r = d$J,
  nI_r = d$I,
  nK_r = d$K,
  xmu_r = d$priors$x$mu,
  xsigma_r = d$priors$x$sigma,
  x_init_r = d$x_init,
  betamu_r = d$priors$beta$mu,
  betasigma_r = d$priors$beta$sigma,
  beta_init_r = d$beta_init,
  verbose_r = TRUE,
  maxit_r = maxit,
  thresh_r = thresh,
  checkfreq_r = checkfreq,
  threads_r = 1,
  p_r = d$token_count,
  prop_ns_r = prop_ns,
  cx_0_r = cx_0,
  dx_0_r = dx_0,
  cb_0_r = cb_0,
  db_0_r = db_0)
  
  x <- bwe_out$x
  colnames(x) <- d$mat_Y@Dimnames[[1]]
  
  beta <- bwe_out$beta
  colnames(beta) <- d$mat_Y@Dimnames[[2]]

return(list(data_prep = d, 
  x = x, 
  beta = beta, 
  var_x = bwe_out$var_x, 
  var_beta = bwe_out$var_beta,
  xsigma = bwe_out$xsigma,
  betasigma = bwe_out$betasigma,
  cx = bwe_out$cx, dx = bwe_out$dx,  
  cb = bwe_out$cb, db = bwe_out$db, 
  convergence = bwe_out$convergence,
  elbo = bwe_out$elbo,
  seed_used = s))
}
