#' Scales Word Embeddings
#' 
#' Scales the word embeddings produced from \code{bwe} according to either one reference word, or a pair of opposite words.
#' 
#' @param target_word Character vector. A single word, used to set the scale the embeddings. If two words are provided, it is assumed they are opposites, they anchor each end of the embedding scaling.
#' @param bwe_object List. Output from \code{fit_bwe}.
#' @param trim Logical. Should common stopwords be removed from the list of scaled words? Setting this to \code{TRUE} improves performance and interpretability. Defaults to \code{TRUE}. 
#' @import data.table
#' @import Matrix
#' @import MASS
#' @import quanteda
#' @return A \code{data.table}. Contains each word, and its position relative to the scaling word(s).

scaling_word_embeddings <- function(target_word, bwe_object, trim = TRUE){
  if(length(target_word) > 2){
    warning("Length of target_word greater than two. Only using first two words.")
  }
  
  x <- bwe_object$x
  
  if(trim == TRUE){
    stops <- quanteda::stopwords("en")
    cols_to_keep <- !(colnames(x) %in% stops)
    x <- x[, cols_to_keep]
  }
  
  
  if(length(target_word) == 1){
    scaled_embedding <- single_anchor(x_matrix = x, target = target_word)
  }
  
  if(length(target_word) == 2){
    scaled_embedding <- double_anchor(
      x_matrix = x, 
      target_word1 = target_word[1],
      target_word2 = target_word[2])
  }
  
  return(scaled_embedding)
}

single_anchor <- function(x_matrix, target){
  x_matrix <- t(x_matrix)
  
  x_vec <- x_matrix[ which(row.names(x_matrix) == target), ]
  word_vec <- matrix(c(x_vec, 1), nrow = 1)
  target_vec <- matrix(0, nrow = ncol(word_vec), ncol = 1)
  target_vec[1, 1] <- 1
  inv_A <- ginv(word_vec)
  foo <- inv_A %*% t(target_vec)
  x2 <- cbind(x_matrix, 1)
  x3 <- x2 %*% foo
  word_dim <- data.table(words = row.names(x_matrix), dim = x3[, 1])
  return(word_dim)
}


double_anchor <- function(x_matrix, target_word1, target_word2){
  x_matrix <- t(x_matrix)
  
  x_vec1 <- x_matrix[ which(row.names(x_matrix) == target_word1), ]
  x_vec2 <- x_matrix[ which(row.names(x_matrix) == target_word2), ]
  x_vec <- rbind(x_vec1, x_vec2)
  word_vec <- cbind(x_vec, 1)
  target_vec <- matrix(0, nrow = ncol(word_vec), ncol = 2)
  target_vec[1, 1] <- 1
  target_vec[1, 2] <- -1
  inv_A <- ginv(word_vec)
  foo <- inv_A %*% t(target_vec)
  x2 <- cbind(x_matrix, 1)
  x3 <- x2 %*% foo
  word_dim <- data.table(words = row.names(x_matrix), dim = x3[, 1])
  return(word_dim)
}
