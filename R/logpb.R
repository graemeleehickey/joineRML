#' @keywords internal
#' @importFrom mvtnorm dmvnorm
logpb <- function(b, theta, data) {

  r <- data$r
  q <- data$q
  K <- data$K
  nk <- data$nk
  y <- data$y
  X <- data$X
  Z <- data$Z
  v <- data$v
  Z.fail <- data$Z.fail
  IW.fail <- data$IW.fail
  tj.ind <- data$tj.ind

  beta <- theta$beta
  gamma <- theta$gamma
  sigma2 <- theta$sigma2
  haz <- theta$haz
  D <- theta$D
  if (sum(r) == 1) {
    D <- as.matrix(D)
  }

  if (length(b) != sum(r)) {
    stop("Incorrect length of b")
  }

  # Expanded gamma_y (repeated for each random effect term)
  if (q > 0) {
    gamma.scale <- diag(rep(gamma[-(1:q)], r), ncol = sum(r))
  } else {
    gamma.scale <- diag(rep(gamma, r), ncol = sum(r))
  }

  # f(b | theta)
  pb <- mvtnorm::dmvnorm(b, mean = rep(0, sum(r)), sigma = D, log = TRUE)

  # f(y | b, theta)
  XbetaZb <- as.vector((X %*% beta) + (Z %*% b))
  Sigma <- diag(x = rep(sigma2, nk), ncol = sum(nk))
  py.b <- mvtnorm::dmvnorm(y, mean = XbetaZb, sigma = Sigma, log = TRUE)

  # f(T > t | b, theta)
  IZ <- t(IW.fail %*% Z.fail)
  W2 <- t(b) %*% gamma.scale %*% IZ
  if (q > 0) {
    W2 <- W2 + as.numeric(t(v) %*% gamma[1:q])
  }
  W2 <- as.vector(W2)
  if (tj.ind > 0) {
    pt.b <- -sum(haz[1:tj.ind] * exp(W2))
  } else {
    pt.b <- 0 # obs times before failure time => survival prob = 1
  }

  out <- pt.b + py.b + pb
  return(out)

}


#' @keywords internal
b_mode <- function(theta, data) {

  out <- optim(par = rep(0, sum(data$r)),
               fn = logpb,
               theta = theta,
               data = data,
               control = list(fnscale = -1),
               method = "BFGS",
               hessian = TRUE)

  return(out)

}


#' @keywords internal
b_metropolis <- function(theta.samp, delta.prop, sigma.prop, b.curr, data.t) {

  accept <- 0

  # Draw b from proposal distribution
  b.prop <- mvtnorm::rmvt(n = 1,
                          delta = delta.prop,
                          sigma = sigma.prop,
                          df = 4)
  b.prop <- as.vector(b.prop)

  # Metropolis-Hastings acceptance
  log.a1 <- logpb(b.prop, theta.samp, data.t) - logpb(b.curr, theta.samp, data.t)

  dens.curr <- mvtnorm::dmvt(x = b.curr,
                             delta = delta.prop,
                             sigma = sigma.prop,
                             df = 4,
                             log = TRUE)
  dens.prop <- mvtnorm::dmvt(x = b.prop,
                             delta = delta.prop,
                             sigma = sigma.prop,
                             df = 4,
                             log = TRUE)
  log.a2 <- dens.curr - dens.prop
  a <- min(exp(log.a1 - log.a2), 1)
  randu <- runif(1)
  if (randu <= a) {
    b.curr <- b.prop
    accept <- 1
  }


  out <- list("b.curr" = b.curr, "accept" = accept)
  return(out)

}

