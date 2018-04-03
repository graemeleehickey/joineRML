#' @keywords internal
#' @importFrom randtoolbox sobol halton
#' @importFrom MASS mvrnorm
stepEM <- function(theta, l, t, z, nMC, verbose, gammaOpt, pfs, type) {

  # Input parameter estimates
  D <- theta$D
  beta <- theta$beta
  sigma2 <- theta$sigma2
  haz <- theta$haz
  gamma <- theta$gamma

  # Multivariate longitudinal data
  yi <- l$yi
  Xi <- l$Xi
  XtX.inv <- l$XtX.inv
  Xtyi <- l$Xtyi
  XtZi <- l$XtZi
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

  # MVN mean vector for [y | b]
  Mi <- mapply(function(a, z, s, y, X) {
    as.vector(a %*% (z %*% s %*% (y - X %*% beta)))
  },
  a = Ai, z = Zit, s = Sigmai.inv, y = yi, X = Xi,
  SIMPLIFY = FALSE)

  # Monte Carlo sample of [b | y]
  if (type == "montecarlo") {
    bi.y <- mapply(function(m, a) {
      MASS::mvrnorm(nMC, mu = m, Sigma = a)
    },
    m = Mi, a = Ai,
    SIMPLIFY = FALSE)
  } else if (type == "antithetic") {
    bi.y <- bSim(floor(nMC / 2),  Mi, Ai)
  } else if (type %in% c("sobol", "halton")) {
    if (type == "sobol") {
      Zq <- randtoolbox::sobol(nMC, dim = sum(r), normal = TRUE, scrambling = 1)
    } else if (type == "halton") {
      Zq <- randtoolbox::halton(nMC, dim = sum(r), normal = TRUE)
    }
    bi.y <- mapply(function(m, a) {
      C <- chol(a)
      matrix(rep(m, nMC), nrow = nMC, byrow = TRUE) + (Zq %*% C)
    },
    m = Mi, a = Ai,
    SIMPLIFY = FALSE)
  }
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
    gamma.scale <- diag(rep(gamma[-(1:q)], r), ncol = sum(r))
  } else {
    gamma.scale <- diag(rep(gamma, r), ncol = sum(r))
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
  rm(expW) # don't need this anymore (large memory object)

  # lambda0(t) for profile score function of beta
  haz.hat <- hazHat(expvstargam, pb.yt, nev)
  haz.hat <- as.vector(haz.hat)

  if (gammaOpt == "GN") {
    gDelta <- gammaUpdate_approx(bi.y, Zit.fail, expvstargam, pb.yt, haz.hat,
                                 V, survdat2.list, K, q, nev.uniq)$gDelta
  } else {
    gDelta <- gammaUpdate(bi.y, Zit.fail, expvstargam, pb.yt, haz.hat,
                          V, survdat2.list, K, q, nev.uniq, nev)$gDelta
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
  rr <- mapply(function(x1, x2, b) {
    x1 - (x2 %*% b)
  },
  x1 = Xtyi, x2 = XtZi, b = Eb)
  rr.sum <- rowSums(rr)

  beta.new <- as.vector(XtX.inv %*% rr.sum)
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
      t(residFixed) %*% (residFixed - 2*(z %*% b.k)) + sum(diag(crossprod(z) %*% bbT.k))
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
    gamma.new.scale <- diag(rep(gamma.new[-(1:q)], r), ncol = sum(r))
  } else {
    gamma.new.scale <- diag(rep(gamma.new, r), ncol = sum(r))
  }

  haz.new <- lambdaUpdate(bi.y, IW.fail, Zi.fail, pb.yt, V,
                          gamma.new.scale, gamma.new, q, nev, survdat2.list)
  haz.new <- as.vector(haz.new)

  theta.new <- list("D" = D.new, "beta" = beta.new, "sigma2" = sigma2.new,
                    "haz" = haz.new, "gamma" = gamma.new)

  t3 <- Sys.time()

  if (verbose && !pfs) {
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

  ## These are only calculated when pfs = TRUE

  if (pfs) {

    # Posterior means + variances of random effects

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

    out$Eb <- Eb.flat
    out$Vb <- Vb

    #-----------------------------------------------------

    # Approximate standard errors

    m <- list()
    m$Sigmai.inv <- Sigmai.inv
    m$Eb <- Eb
    m$EbbT <- EbbT
    m$bi.y <- bi.y
    m$expvstargam <- expvstargam
    m$pb.yt <- pb.yt
    m$haz.hat <- haz.hat

    out$Hessian <- hessian(theta, l, t, z, m)

  }

  return(out)

}
