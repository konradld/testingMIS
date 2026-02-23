#' Frisch-Waugh-Lovell projection: residualizes y and X1 on X2
fwl <- function(y, X1, X2) {
  Q <- qr(as.matrix(X2))
  cbind(Yp = y - qr.fitted(Q, y), Xp = as.matrix(X1) - qr.fitted(Q, as.matrix(X1)))
}