#' Internal function for performing a single iteration of the MCEM algorithm
#'
#' Also takes logical arguments \code{ll} and \code{se.approx}, which calculates
#' the log-likelihood (and posterior mean and variance of the random effects)
#' and approximate standard errors, respectively. When either of these arguments
#' are \code{TRUE}, the function does not return the maaximizer from the MCEM
#' algorithm iteration, and instead reports the called for post-fit statistics.
#'
#' @keywords internal
stepEM <- function(theta, l, t, z, nMC, verbose, gammaOpt, postRE, se.approx) {

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

  # Inverse-Sigma_i (error precision matrix; diagonal matrix)
  Sigmai.inv <- lapply(nik, function(i) {
    diag(x = rep(1 / sigma2, i), ncol = sum(i))
  })

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
  #     a <- fast_nearPD(a)
  #   } else {
  #     a
  #   }
  # })

  # MVN mean vector for [y | b]
  Mi <- mapply(function(a, z, s, y, X) {
    as.vector(a %*% (z %*% s %*% (y - X %*% beta)))
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
  if (sum(r) > 1) {
    if (q > 0) {
      gamma.scale <- diag(rep(gamma[-(1:q)], r))
    } else {
      gamma.scale <- diag(rep(gamma, r))
    }
  } else {
    gamma.scale <- gamma[length(gamma)] # just a single gamma_y
  }

  # exp{W(tj, b)}
  IZ <- mapply(function(x, y) {
    t(x %*% y)
  },
  x = IW.fail, y = Zi.fail,
  SIMPLIFY = FALSE)
  # subjects who are censored before first failure time do not contribute anything
  # -> this information is captured through expW
  expW <- expWArma(IZ, bi.y, gamma.scale, survdat2.list)

  # log{f(T, delta | b)}
  logfti <- mapply(function(w, v, h) {
    H <- as.vector(w %*% haz[1:ncol(w)]) * exp(v) # cummulative hazard
    if (h$delta == 1) { # event
      (log(haz[ncol(w)]) + v + log(w[, ncol(w)])) - H
    } else { # non-event
      -H
    }
  },
  w = expW, v = Vtgamma, h = survdat2.list,
  SIMPLIFY = FALSE)
  fti <- lapply(logfti, exp) # f(T, delta | b)

  # Expectation denominator
  den <- lapply(fti, mean)

  # f(T, delta | b) / den
  pb.yt <- mapply(function(f, d) {
    f / d
  },
  f = fti, d = den,
  SIMPLIFY = FALSE)

  t1 <- Sys.time()

  #*****************************************************
  # E-step starts here
  #*****************************************************

  # E[b]
  Eb <- mapply(function(b, pb) {
    colMeans(b * pb)
  },
  b = bi.y, pb = pb.yt,
  SIMPLIFY = FALSE)

  # E[bb^T]
  EbbT <- mapply(function(b, pb) {
    crossprod(b, (b * pb)) / nrow(b)
  },
  b = bi.y, pb = pb.yt,
  SIMPLIFY = FALSE)

  # exp{v %*% gamma_v + W(tj)}
  expvstargam <- mapply(function(w, v) {
    w * exp(v)
  },
  w = expW, v = Vtgamma,
  SIMPLIFY = FALSE)

  # lambda0(t) for profile score function of beta
  haz.hat <- hazHat(expvstargam, pb.yt, nev)
  haz.hat <- as.vector(haz.hat)

  if (gammaOpt == "GN") {
    gDelta <- gammaUpdate_approx(bi.y, Zit.fail, expvstargam, pb.yt, haz.hat,
                                 V, survdat2.list, K, q, nev.uniq)$gDelta
  } else {
    gDelta <- gammaUpdate(bi.y, Zit.fail, expvstargam, pb.yt, haz.hat,
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
  XtX <- mapply(function(xt, x) {
    xt %*% x
  },
  xt = Xit, x = Xi,
  SIMPLIFY = FALSE)
  XtX.sum <- Reduce("+", XtX)

  rr <- mapply(function(xt, y, z, b) {
    xt %*% (y - (z %*% b))
  },
  xt = Xit, y = yi, z = Zi, b = Eb)
  rr.sum <- rowSums(rr)

  beta.new <- solve(XtX.sum, rr.sum)
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

  haz.new <- lambdaUpdate(bi.y, IW.fail, Zi.fail, pb.yt, V,
                          gamma.new.scale, gamma.new, q, nev, survdat2.list)
  haz.new <- as.vector(haz.new)

  theta.new <- list("D" = D.new, "beta" = beta.new, "sigma2" = sigma2.new,
                    "haz" = haz.new, "gamma" = gamma.new)

  t3 <- Sys.time()

  if (verbose && !postRE) {
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

  #*********************************************************
  # Post model fit processing: log-likehood + posterior REs
  #*********************************************************

  ## Observed log-likelihood (for iteration t-1)

  fy <- mapply(function(y, x, z, zt, s, nik) {
    # marginal (observed data) likelihood for long. data
    r <- y - (x %*% beta)
    v <- s + z %*% D %*% zt
    vinv <- solve(v)
    -0.5 * (sum(nik) * log(2*pi) +
              as.numeric(determinant(v, logarithm = TRUE)$modulus) +
              colSums(t(r) %*% vinv %*% r))

  },
  y = yi, x = Xi, z = Zi, zt = Zit, s = Sigmai, nik = nik)

  ll <- sum(fy + log(unlist(den)))
  out <- list("theta.new" = theta.new, "ll" = ll)

  #-----------------------------------------------------

  ## Posterior means + variances of random effects

  if (postRE) {

    # Var(b)
    Vb <- mapply(function(b, pb, mu) {
      v <- (crossprod(b, (b * pb)) / nrow(b)) - tcrossprod(mu)
      rownames(v) <- colnames(v) <- colnames(D)
      v
    },
    b = bi.y, pb = pb.yt, mu = Eb,
    SIMPLIFY = "array")

    # E[b]
    Eb.flat <- simplify2array(Eb)
    if (sum(r) > 1) {
      Eb.flat <- t(Eb.flat)
    } else {
      Eb.flat <- as.matrix(Eb.flat, ncol = 1)
    }
    colnames(Eb.flat) <- colnames(D)

    out$Eb = Eb.flat
    out$Vb = Vb

  }

  #*********************************************************
  # Approximate standard errors
  #*********************************************************

  if (se.approx) {

    m <- list()
    m$Sigmai.inv <- Sigmai.inv
    m$Eb <- Eb
    m$EbbT <- EbbT
    m$bi.y <- bi.y
    m$expvstargam <- expvstargam
    m$pb.yt <- pb.yt
    m$haz.hat <- haz.hat

    out$ses <- approxSE(theta, l, t, z, m)

  }

  return(out)

}
