#' Frisch-Waugh-Lovell projection: residualizes y and X1 on X2
#'
#' @noRd
fwl <- function(y, X, Z) {
  Q <- qr(as.matrix(Z))
  cbind(
    Yp = y - qr.fitted(Q, y),
    Xp = as.matrix(X) - qr.fitted(Q, as.matrix(X))
  )
}

#' Splits a vector into equal-sized blocks (drops remainder)
#'
#' @noRd
make_blocks <- function(X, block_size) {
  n <- length(X) %/% block_size
  matrix(X[seq_len(n * block_size)], nrow = block_size, ncol = n)
}
