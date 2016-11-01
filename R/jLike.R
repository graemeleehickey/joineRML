#' Calculate the log-likelihood and posterior random effects at the maximizer
#'
#' @keywords internal
jLike <- function(theta, l, t, z, nMC) {

  # MLE parameter estimates from EM algorithm
  D <- theta$D
  beta <- theta$beta
  sigma2 <- theta$sigma2
  haz <- theta$haz
  gamma <- theta$gamma

  # Multivariate longitudinal data
  yi <- l$yi
  Xi <- l$Xi
  Zi <- l$Zi
  Zit <- l$Zit
  nik <- l$nik
  r <- l$r # vector of random effect dims

  # Time-to-event data
  V <- t$V
  survdat2.list <- t$survdat2.list
  q <- t$q

  # Covariate data for W(u, b)
  Zi.fail <- z$Zi.fail
  IW.fail <- z$IW.fail

  #*****************************************************
  # Monte Carlo set-up
  #*****************************************************

  # Sigma_i (error covariance matrix; diagonal matrix)
  Sigmai <- lapply(nik, function(i) {
    diag(x = rep(sigma2, i), ncol = sum(i))
  })
  Sigmai.inv <- lapply(Sigmai, solve)

  # MVN covariance matrix for [b | y]
  Dinv <- solve(D)
  Ai <- mapply(FUN = function(zt, s, z) {
    solve((zt %*% s %*% z) + Dinv)
  },
  z = Zi, zt = Zit, s = Sigmai.inv,
  SIMPLIFY = FALSE)

  # MVN mean vector for [y | b]
  Mi <- mapply(function(a, z, s, y, X) {
    a %*% (z %*% s %*% (y - X %*% beta))
  },
  a = Ai, z = Zit, s = Sigmai.inv, y = yi, X = Xi,
  SIMPLIFY = FALSE)

  # Monte Carlo sample of [b | y]
  bi.y <- bSim(floor(nMC / 2),  Mi, Ai)
  names(bi.y) <- names(Ai)

  # Calculation of t(v) %*% gamma_v
  if (q > 0) {
    Vtgamma <- mapply(function(v) {
      as.numeric(v %*% gamma[1:q])
    },
    v = V,
    SIMPLIFY = FALSE)
  } else {
    Vtgamma <- V
  }

  # Expanded gamma_y (repeated for each random effect term)
  if (q > 0) {
    gamma.scale <- diag(rep(gamma[-(1:q)], r))
  } else {
    gamma.scale <- diag(rep(gamma, r))
  }

  # exp{W(tj, b)}
  IZ <- mapply(function(x, y) t(x %*% y),
               x = IW.fail, y = Zi.fail,
               SIMPLIFY = FALSE)
  gamma.b <- lapply(bi.y, function(x) x %*% gamma.scale)
  expW <- mapply(function(x, y, h) {
    exp(y %*% x) * ifelse(h$tj.ind > 0, 1, 0) # subjects who are censored before first
  },                                          # failure time do not contribute anything
  x = IZ, y = gamma.b, h = survdat2.list)

  # Cummulative hazard function
  Hi <- mapply(function(w, v) {
    as.vector(w %*% haz[1:ncol(w)]) * exp(v)
  },
  w = expW, v = Vtgamma,
  SIMPLIFY = FALSE)

  # f(T, delta | b)
  fti <- mapply(function(h, H, v, w) {
    (haz[ncol(w)] * exp(v) * w[, ncol(w)])^(h$delta) * exp(-H)
  },
  h = survdat2.list, H = Hi, w = expW, v = Vtgamma,
  SIMPLIFY = FALSE)

  # Expectation denominator
  den <- lapply(fti, mean)

  #*****************************************************
  # Observed log-likelihood
  #*****************************************************

  # f(y): longitudinal data marginal (observed data) likelihood
  fy <- mapply(function(y, x, z, zt, s, nik) {
    r <- y - (x %*% beta)
    v <- s + z %*% D %*% zt
    vinv <- solve(v)
    -0.5 * (sum(nik) * log(2*pi) +
              as.numeric(determinant(v, logarithm = TRUE)$modulus) +
              colSums(t(r) %*% vinv %*% r))

  },
  y = yi, x = Xi, z = Zi, zt = Zit, s = Sigmai, nik = nik)

  ll <- sum(fy + log(unlist(den)))

  #*****************************************************
  # Posterior means of random effects
  #*****************************************************

  # E[b]
  Eb <- mapply(function(b, f, d) {
    colMeans(b * f) / d
  },
  b = bi.y, f = fti, d = den,
  SIMPLIFY = FALSE)

  # Var(b)
  Vb <- mapply(function(b, f, d, mu) {
    v <- crossprod(b, (b * f)) / (nrow(b) * d) - tcrossprod(mu)
    rownames(v) <- colnames(v) <- colnames(D)
    v
  },
  b = bi.y, f = fti, d = den, mu = Eb,
  SIMPLIFY = "array")

  Eb <- t(simplify2array(Eb))
  colnames(Eb) <- colnames(D)

  return(list(ll = ll, Eb = Eb, Vb = Vb))

}
