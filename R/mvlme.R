#' @keywords internal
mvlme <- function(thetaLong, l, z, tol.em, verbose) {

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
  n <- l$n    # number of subjects
  p <- l$p    # vector of fixed effect dims
  r <- l$r    # vector of random effect dims
  K <- l$K    # number of longitudinal markers
  nk <- l$nk  # vector of number of observations per outcome

  delta <- 1
  while (delta > tol.em) {

    # Input parameter estimates
    D <- thetaLong$D
    beta <- thetaLong$beta
    sigma2 <- thetaLong$sigma2

    #*****************************************************
    # E-step
    #*****************************************************

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
    Eb <- mapply(function(a, z, s, y, X) {
      as.vector(a %*% (z %*% s %*% (y - X %*% beta)))
    },
    a = Ai, z = Zit, s = Sigmai.inv, y = yi, X = Xi,
    SIMPLIFY = FALSE)

    # E[bb^T]
    EbbT <- mapply(function(v, e) {
      v + tcrossprod(e)
    },
    v = Ai, e = Eb,
    SIMPLIFY = FALSE)

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

    thetaLong.new <- list("D" = D.new, "beta" = beta.new, "sigma2" = sigma2.new)

    # Relative parameter change
    delta <- sapply(c("D", "beta", "sigma2"), function(i) {
      abs(thetaLong[[i]] - thetaLong.new[[i]]) / (abs(thetaLong[[i]]) + 1e-03)
    })
    delta <- max(unlist(delta))

    thetaLong <- thetaLong.new

    if (verbose) {
      print(thetaLong.new)
    }

  }

  return(thetaLong.new)

}
