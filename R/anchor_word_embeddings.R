#' Anchor Word Embeddings
#' 
#' Anchors and identifies the embeddings produced from \code{bwe}.
#' 
#' @param anchoring_word Character. A single word, used to set the anchorings.
#' @param bwe_object List. Output from \code{fit_bwe}.
#' @param similarity Character. Should similarity be measured only by cosine similiarity, or should similarity be weighted by the log frequency of words.
#' @param trim Boolean. Should common stopwords be removed from the list of possible anchors? Setting this to \code{TRUE} improves performance and interpretability. Defaults to \code{TRUE}. 
#' @import data.table
#' @import Matrix
#' @import quanteda
#' @import stringr
#' @return A list. Contains a matrix of the anchored embeddings, and vector of the words used as anchors.
#' @export

anchor_word_embeddings <- function(anchoring_word, bwe_object, 
  similarity = c("cosine", "weighted cosine"), trim = TRUE){
  data_to_use <- bwe_object$data_prep
  words <- data_to_use$mat_Y@Dimnames[[1]]
  word_frq <- data_to_use$token_count_int
  word_frq <- word_frq[ order(word_frq$N, decreasing = TRUE)]
  word_frq[, log_freq := log(N)]
  # word_frq <- word_frq[1:floor((1-drop) * nrow(word_frq))]
  if(!(anchoring_word[1] %in%word_frq$words )){
    stop(paste(anchoring_word[1], "doesn't appear"))
  }
  nK <- nrow(bwe_object$x)
  word_index_1 <- grep(pattern = paste0("^", anchoring_word[1], "$"), x = words)
  
  mat1 <- bwe_object$x
  mat1 <- t(bwe_object$x)
  mat_magnitudes <- rowSums(mat1^2)
  
  anchoring_word_indices <- matrix(NA_integer_, nrow = nK + 1, ncol = 1)
  anchoring_words <- matrix(NA_character_,  nrow = nK + 1, ncol = 1)
  anchoring_word_indices[1, ] <- word_index_1
  anchoring_words[1, ] <- anchoring_word[1]
  
  # given initial word, return words which are opposites in cosine space
  for(i in 1:(nK+1)){
    vec1 <- bwe_object$x[, anchoring_word_indices[i] ]
    vec1 <- Matrix(vec1, nrow = 1, ncol = length(vec1))
    vec_magnitudes <- rowSums(vec1^2)
    sim <- matrix(t(tcrossprod(vec1, mat1)/
        (sqrt(tcrossprod(vec_magnitudes, mat_magnitudes)))))
    sim2 <- data.table(words = words, cosine_sim = sim)
    sim2 <- sim2[ words %in% word_frq$words ]
    sim2 <- merge(sim2, word_frq, by = "words")
    sim2[, weight_cos := (cosine_sim.V1) * log_freq]
    if(trim == TRUE){
      sim2 <- sim2[ nchar(words) > 3 ]
      sim2 <- sim2[ !(words %in% quanteda::stopwords("en"))]
      if(similarity == "weighted cosine"){
        w <- sim2[ order(weight_cos) ]  
        w <- w[ !(words %in% anchoring_words) ]
      } else{
        w <- sim2[ order(cosine_sim.V1) ]  
        w <- w[ !(words %in% anchoring_words) ]
      }
    }else{
      if(similarity == "weighted cosine"){
        w <- sim2[ order(weight_cos) ]  
        w <- w[ !(words %in% anchoring_words) ]
      } else{
        w <- sim2[ order(cosine_sim.V1) ]  
        w <- w[ !(words %in% anchoring_words) ]
      }
    }
    
    if(i != (nK+1)){
      if(is.na(anchoring_word[i + 1])){
        anchoring_words[i + 1, ] <- w[1, ]$words
      } else{
        anchoring_words[i + 1, ] <- anchoring_word[i + 1]
      }
      word_regex <- paste0("^", anchoring_words[i + 1, ], "$")
      possible_matches <- words[grep(pattern = word_regex, x = words)]
      word_idx <- grep(
        pattern = word_regex, x = words)[
          anchoring_words[i + 1, ] == possible_matches]
      anchoring_word_indices[i+1, ] <- word_idx
    }
  }
  
  # creating the anchor points from the words of choice
  x_anchors <- bwe_object$x[, anchoring_word_indices]
  positive_anchors <- seq(1, nK+1, 2)
  negative_anchors <- seq(2, nK+1, 2)
  if(length(anchoring_word_indices) %% 2 != 0){
    positive_anchors <- positive_anchors[1:(length(positive_anchors)-1)]
  }
  if(length(anchoring_word_indices) %% 2 == 0){
    negative_anchors <- negative_anchors[1:(length(negative_anchors)-1)]
  }

  
  
  target_vectors <- x_anchors * 0
  for(i in 1:length(positive_anchors)){
    target_vectors[positive_anchors[i], positive_anchors[i] ] <- 1
  }
  for(i in 1:length(negative_anchors)){
    target_vectors[ negative_anchors[i], negative_anchors[i] ] <- -1
  }
  word_embeddings <- t(bwe_object$x)
  x_anchors <- t(x_anchors)
  target_vectors <- t(target_vectors)
  # affine transformation (from pscl:::affineTrans )
  x_anchors_0 <- cbind(x_anchors, 1)
  # zero_mat <- 0 * x_anchors_0
  # we need a block diagonal matrix here, repeat anchors nK times
  A <- bdiag(rep(list(x_anchors_0), nK)) 
  b <- matrix(target_vectors, ncol = 1)
  foo <- solve(A, tol = 1e-18) %*% b
  foo <- matrix(foo, nrow = nK + 1)
  word_embeddings <- cbind(word_embeddings, 1)
  identified_words <- word_embeddings %*% foo
  identified_words <- t(identified_words)
  colnames(identified_words) <- words
  return(list(identified_words = identified_words, 
    anchoring_words = anchoring_words))
}