#' Internal function for performing a single iteration of the MCEM algorithm
#'
#' @keywords internal
stepEM <- function(theta, l, t, z, nMC, verbose, approxInfo) {

  # Input parameter estimates
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

  t0 <- Sys.time()

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

  # # If Ai is non-PSD, approximate Ai by nearest PD matrix
  # Ai <- lapply(Ai, FUN = function(a) {
  #   if (any(eigen(a, symmetric = TRUE)$values < 0)) {
  #     print("Non-PSD matrix detected!")
  #     a <- nearPD(a)
  #   } else {
  #     a
  #   }
  # })

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

  t1 <- Sys.time()

  #*****************************************************
  # E-step starts here
  #*****************************************************

  # E[b]
  Eb <- mapply(function(b, f, d) {
    colMeans(b * f) / d
  },
  b = bi.y, f = fti, d = den,
  SIMPLIFY = FALSE)

  # E[bb^T]
  EbbT <- mapply(function(b, f, d) {
    crossprod(b, (b * f)) / (nMC * d)
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

  if (approxInfo) {
    gDelta <- gammaUpdate_approx(bi.y, Zit.fail, expvstargam, fti, den, haz.hat,
                                 V, survdat2.list, K, q, nev.uniq)$gDelta
  } else {
    gDelta <- gammaUpdate(bi.y, Zit.fail, expvstargam, fti, den, haz.hat,
                          V, survdat2.list, K, q, nev.uniq)$gDelta
  }

  t2 <- Sys.time()

  #*****************************************************
  # M-step starts here
  #*****************************************************

  # D
  D.new <- Reduce("+", EbbT) / n
  rownames(D.new) <- colnames(D.new) <- rownames(D)

  #-----------------------------------------------------

  # beta
  XtSX <- mapply(function(xt, x) {
    xt %*% x
  },
  xt = Xit, x = Xi,
  SIMPLIFY = FALSE)
  XtSX.sum <- Reduce("+", XtSX)

  rr <- mapply(function(xt, y, z, b) {
    xt %*% (y - z %*% b)
  },
  xt = Xit, y = yi, z = Zi, b = Eb)
  rr.sum <- rowSums(rr)

  beta.new <- solve(XtSX.sum, rr.sum)
  names(beta.new) <- names(beta)

  #-----------------------------------------------------

  # sigma_k^2
  beta.inds <- cumsum(c(0, p))
  b.inds <- cumsum(c(0, r))
  sigma2.new <- vector(length = K)

  for (k in 1:K) {
    beta.k <- beta.new[(beta.inds[k] + 1):(beta.inds[k + 1])]
    SSq <- mapply(function(y, x, z, b, b2) {
      b.k <- b[(b.inds[k] + 1):(b.inds[k + 1])]
      bbT.k <- b2[(b.inds[k] + 1):(b.inds[k + 1]), (b.inds[k] + 1):(b.inds[k + 1])]
      residFixed <- (y - x %*% beta.k)
      t(residFixed) %*% (residFixed - 2*(z %*% b.k)) + sum(diag((t(z) %*% z) %*% bbT.k))
    },
    y = yik[[k]], x = Xik.list[[k]], z = Zik.list[[k]], b = Eb, b2 = EbbT)
    sigma2.new[k] <- sum(SSq) / nk[[k]]
  }

  names(sigma2.new) <- paste0("sigma2_", 1:K)

  #-----------------------------------------------------

  # gamma
  gamma.new <- gamma + as.vector(gDelta)

  #-----------------------------------------------------

  # lambda0(tj)

  # Expanded gamma_y (repeated for each random effect term)
  # - using the latest EM iteration estimate
  if (q > 0) {
    gamma.new.scale <- diag(rep(gamma.new[-(1:q)], r))
  } else {
    gamma.new.scale <- diag(rep(gamma.new, r))
  }

  haz.new <- lambdaUpdate(bi.y, IW.fail, Zi.fail, fti, V, den,
                          gamma.new.scale, gamma.new, q, nev, survdat2.list)
  haz.new <- as.vector(haz.new)

  #-----------------------------------------------------

  theta.new <- list("D" = D.new, "beta" = beta.new, "sigma2" = sigma2.new,
                    "haz" = haz.new, "gamma" = gamma.new)
  theta = theta.new

  t3 <- Sys.time()

  if (verbose) {
    tdiff1 <- t1 - t0
    cat(paste("Step 1: Time to setup Monte Carlo expectations", round(tdiff1, 2),
              attr(tdiff1, "units"), "\n"))
    tdiff2 <- t2 - t1
    cat(paste("Step 2: Time to perform E-step", round(tdiff2, 2),
              attr(tdiff2, "units"), "\n"))
    tdiff3 <- t3 - t2
    cat(paste("Step 3: Time to perform M-step", round(tdiff3, 2),
              attr(tdiff3, "units"), "\n"))
    tdiff4 <- t3 - t0
    cat(paste("Total time for EM algorithm", round(tdiff4, 2),
              attr(tdiff4, "units"), "\n"))
  }

  return(theta.new)

}
