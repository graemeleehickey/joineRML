#' Approximate standard errors using the empirical information matrix
#'
#' @keywords internal
approxSE <- function(theta, l, t, z, nMC) {

  # MLE parameter estimates from EM algorithm
  D <- theta$D
  beta <- theta$beta
  sigma2 <- theta$sigma2
  haz <- theta$haz
  gamma <- theta$gamma

  # Multivariate longitudinal data
  yi <- l$yi
  Xi <- l$Xi
  Xit <- l$Xit
  Zi <- l$Zi
  Zit <- l$Zit
  nik <- l$nik
  yik <- l$yik
  Xik.list <- l$Xik.list
  Zik.list <- l$Zik.list
  p <- l$p    # vector of fixed effect dims
  r <- l$r    # vector of random effect dims
  K <- l$K    # number of longitudinal markers
  n <- l$n    # number of subjects
  nk <- l$nk  # vector of number of observations per outcome

  # Time-to-event data
  V <- t$V
  survdat2 <- t$survdat2
  survdat2.list <- t$survdat2.list
  q <- t$q
  nev <- t$nev
  nev.uniq <- t$nev.uniq

  # Covariate data for W(u, b)
  Zi.fail <- z$Zi.fail
  Zit.fail <- z$Zit.fail
  Zik.list <- z$Zik.list
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
  # Conditional expectations
  #*****************************************************

  # E[b]
  Eb <- mapply(function(b, f, d) {
    colMeans(b * f) / d
  },
  b = bi.y, f = fti, d = den,
  SIMPLIFY = FALSE)

  # E[bb^T]
  EbbT <- mapply(function(b, f, d) {
    crossprod(b, (b * f)) / (nrow(b) * d)
  },
  b = bi.y, f = fti, d = den,
  SIMPLIFY = FALSE)

  # E[exp{W(tj, b)}]
  EexpW <- EexpWArma(expW, fti, den)
  names(EexpW) <- names(Ai)

  # exp{v %*% gamma_v + W(tj)}
  expvstargam <- mapply(function(w, v) {
    w * exp(v)
  },
  w = expW, v = Vtgamma,
  SIMPLIFY = FALSE)

  # lambda0(t) for profile score function of beta
  haz.hat <- hazHat(expvstargam, fti, den, nev)
  haz.hat <- as.vector(haz.hat)

  #*****************************************************
  # Complete data score components
  #*****************************************************

  # beta

  sbeta <- mapply(function(xt, x, sinv, y, z, b) {
    (xt %*% sinv) %*% (y - x %*% beta - z %*% b)
  },
  xt = Xit, x = Xi, sinv = Sigmai.inv, y = yi, z = Zi, b = Eb,
  SIMPLIFY = TRUE)

  rownames(sbeta) <- names(beta)

  #-----------------------------------------------------

  # Dinv

  sDinv <- lapply(EbbT, function(b2) {
    0.5 * (2*D - diag(D)) - 0.5 * (2 * b2 - diag(b2))
  })
  sDinv <- sapply(sDinv, function(d) d[lower.tri(d, diag = TRUE)])

  ltri <- lower.tri(D, diag = TRUE)
  rownames(sDinv) <- paste0("Dinv_", row = row(D)[ltri], ",", col = col(D)[ltri])

  #-----------------------------------------------------

  # gamma

  sgamma <- gammaUpdate(bi.y, Zit.fail, expvstargam, fti, den, haz.hat,
                        V, survdat2.list, K, q, nev.uniq)$scorei

  rownames(sgamma) <- names(gamma)

  #-----------------------------------------------------

  # sigma2

  beta.inds <- cumsum(c(0, p))
  b.inds <- cumsum(c(0, r))

  ssigma2 <- matrix(nrow = K, ncol = n)
  for (k in 1:K) {
    beta.k <- beta[(beta.inds[k] + 1):(beta.inds[k + 1])]
    ssigma2[k, ] <- mapply(function(y, x, z, b, b2, nik) {
      b.k <- b[(b.inds[k] + 1):(b.inds[k + 1])]
      bbT.k <- b2[(b.inds[k] + 1):(b.inds[k + 1]), (b.inds[k] + 1):(b.inds[k + 1])]
      residFixed <- (y - x %*% beta.k)
      resids <- t(residFixed) %*% (residFixed - 2*(z %*% b.k)) + sum(diag(t(z) %*% z %*% bbT.k))
      (-0.5 * nik[k] / sigma2[k]) + (0.5 * resids / sigma2[k]^2)
    },
    y = yik[[k]], x = Xik.list[[k]], z = Zik.list[[k]], b = Eb, b2 = EbbT, nik = nik)
  }

  rownames(ssigma2) <- paste0("sigma2_", 1:K)

  #-----------------------------------------------------

  si <- rbind(sDinv, sbeta, ssigma2, sgamma)

  ses <- matrix(0, nrow(si), nrow(si))
  for (j in 1:ncol(si)) ses <- ses + si[, j] %*% t(si[, j])
  # Although RHS term = 0 in theory, in practice with MC integration
  # not all terms are vanishingly small, so we add it in
  ses <- ses - (rowSums(si) %*% t(rowSums(si))) / ncol(si)
  rownames(ses) <- colnames(ses) <- rownames(si)

  return(ses)

}
