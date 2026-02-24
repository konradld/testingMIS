#' Frisch-Waugh-Lovell projection: residualizes y and X1 on X2
#'
#' @noRd
fwl <- function(y, X1, X2) {
  Q <- qr(as.matrix(X2))
  cbind(
    Yp = y - qr.fitted(Q, y),
    Xp = as.matrix(X1) - qr.fitted(Q, as.matrix(X1))
  )
}

#' Splits a vector into equal-sized blocks (drops remainder)
#'
#' @noRd
make_blocks <- function(X, block_size) {
  n <- length(X) %/% block_size
  matrix(X[seq_len(n * block_size)], nrow = block_size, ncol = n)
}
