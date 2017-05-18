#' @keywords internal
S <- function(b, theta, data) {

  # NB: landmarking time recorded in 'data'

  r <- data$r
  q <- data$q
  v <- data$v
  Z.fail <- data$Z.fail
  IW.fail <- data$IW.fail
  tj.ind <- data$tj.ind

  gamma <- theta$gamma
  haz <- theta$haz

  if (length(b) != sum(r)) {
    stop("Incorrect length of b")
  }

  # Expanded gamma_y (repeated for each random effect term)
  if (q > 0) {
    gamma.scale <- diag(rep(gamma[-(1:q)], r), ncol = sum(r))
  } else {
    gamma.scale <- diag(rep(gamma, r), ncol = sum(r))
  }

  IZ <- t(IW.fail %*% Z.fail)
  W2 <- t(b) %*% gamma.scale %*% IZ
  if (q > 0) {
    W2 <- W2 + as.numeric(t(v) %*% gamma[1:q])
  }
  W2 <- as.vector(W2)
  if (tj.ind > 0) {
    H <- sum(haz[1:tj.ind] * exp(W2))
  } else {
    H <- 0
  }

  return(exp(-H))

}
