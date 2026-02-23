#' Splits a vector into equal-sized blocks (drops remainder)
make_blocks <- function(X, block_size) {
  n <- length(X) %/% block_size
  matrix(X[seq_len(n * block_size)], nrow = block_size, ncol = n)
}