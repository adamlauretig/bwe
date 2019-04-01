# Makes a set of priors for x and beta 
make_priors <- function (K = 1, sigma_x = 1, sigma_beta = 1) {
    out <- vector(mode = "list")
    out$x <- list(mu = matrix(rep(0.0, K), nrow = K), sigma = sigma_x^2 * 
        diag(K))
    out$beta <- list(mu = matrix(rep(0.0, K), nrow = K), 
      sigma = sigma_beta^2 * diag(K))
    return(out)
}

# initialized values
make_starts <- function (I, J, K, type = c("zeros", "random")){
  if (type == "zeros") {
    starts <- list(beta = {
      matrix(0.0, nrow = K, ncol = J)
    }, x = {
      matrix(0.0, nrow = K, ncol = I)
    })
  }
  else if (type == "random") {
    starts <- list(beta = {
      matrix(rnorm(J * K) * 10,  nrow = K, ncol = J)
    }, x = {
      matrix(rnorm(I * K) * 1, nrow = K, ncol = I)
    })
  }
  else {
    stop("Unknown type.")
  }
  return(starts)
}


data_prep <- function(dat, 
  text_field, 
  window_to_use,
  minimum_word_count,
  K, 
  char_min,
  start_type, 
  rm_nums,
  stops_list){
  
  if(is.data.frame(dat)){
    if(!(data.table::is.data.table(dat))){
      dat <- data.table::as.data.table(dat)
    }
  } else( stop("dat must be a data.frame or data.table"))
  
  txtfile <- unlist(dat[,..text_field ])

  split_strings <- stringr::str_split(txtfile, pattern = "[.] |[?] |[!] ")
  split_strings <- unlist(split_strings)
  split_strings <- split_strings[nchar(split_strings) > 0]
  all_tokens <- tokens(char_tolower(split_strings), remove_punct = TRUE)
  token_dt <- data.table(unlist(all_tokens))
  token_count <- token_dt[,.N, by = V1]
  token_count <- token_count[ order(token_count$N, decreasing = TRUE) ]
  token_count <- token_count[ N >= minimum_word_count ]
  token_count <- token_count[ nchar(V1) > char_min ]
  if(length(stops_list) > 0){
    token_count <- token_count[ !(V1 %in% stops_list)]
  }
  all_tokens <- tokens_keep(tokens(char_tolower(
    split_strings), 
    remove_punct = TRUE, remove_numbers = rm_nums), token_count$V1)
 w <- as.integer(window_to_use)
  cooc_mat <- fcm(all_tokens, 
    context = "window", count = "boolean", window = w, tri = FALSE, 
    ordered = FALSE)

  cooc_mat@x <- ifelse(cooc_mat@x > 0, 1, 0)
  mat_Y <- cooc_mat
  diag(mat_Y) <- 0

  priors <- make_priors(K = K, sigma_x = 1, 
    sigma_beta = 1)
  
  word_order <- data.table(words = cooc_mat@Dimnames$features, 
    index = 1:length( cooc_mat@Dimnames$features))
  tokens_to_use <- merge(word_order, token_count, by.x = "words", by.y = "V1")
  setorder(tokens_to_use, index)
  tokens_to_use[, word_pr := (N/sum(N))^.75 ]
  
  tokens_to_use <- tokens_to_use[ order(N, decreasing = TRUE) ]
  mat_Y <- mat_Y[ tokens_to_use$words,  tokens_to_use$words]

  if(start_type == "svd"){
    svd_Y <- irlba::irlba(A = mat_Y, nv = K, fastpath = TRUE)
    x_init <- t(svd_Y$u)
    beta_init <- t(svd_Y$v)
  } else{
    init_vals <- make_starts(I = nrow(mat_Y), J = ncol(mat_Y), K = K, type = start_type)
    x_init <- init_vals$x    
    beta_init <- init_vals$beta
  }
  
  list(mat_Y = mat_Y,
    token_count = matrix(tokens_to_use$word_pr, ncol = 1), token_count_int = tokens_to_use[,.(words, N)],
    priors = priors, 
    K = K, I = nrow(mat_Y), J = ncol(mat_Y), 
    x_init = x_init, beta_init = beta_init)
}

