#' Return most similar words
#' 
#' Return the most similar words to a given word, as determined by cosine similarity.
#' 
#' @param vec1 Vector. The vector for a word of interest.
#' @param mat1 Matrix. The \code{x} matrix from either \code{fit_bwe} or \code{anchor_word_embeddings}.
#' @param n_similar Integer. Number of similar words to return. Defaults to 10
#' 
#' @return A vector of length \code{n_similar}.
#' @export

nearest_words <- function(vec1, mat1, n_similar = 10){
  word_labels <- colnames(mat1)
  vec1 <- Matrix(vec1, nrow = 1, ncol = length(vec1))
  mat1 <- t(mat1)
  mat_magnitudes <- rowSums(mat1^2)
  vec_magnitudes <- rowSums(vec1^2)
  sim <- (t(tcrossprod(vec1, mat1)/
      (sqrt(tcrossprod(vec_magnitudes, mat_magnitudes)))))
  sim2 <- matrix(sim, dimnames = list(word_labels))
  
  w <- sim2[order(-sim2),,drop = FALSE]
  return(w[1:n_similar,])
}