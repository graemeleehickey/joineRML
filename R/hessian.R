#' @keywords internal
hessian <- function(theta, l, t, z, m) {

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

  # MCEM calculations
  Sigmai.inv <- m$Sigmai.inv
  Eb <- m$Eb
  EbbT <- m$EbbT
  bi.y <- m$bi.y
  expvstargam <- m$expvstargam
  pb.yt <- m$pb.yt
  haz.hat <- m$haz.hat

  #*****************************************************
  # Complete data score components
  #*****************************************************

  # beta

  sbeta <- mapply(function(x, sinv, y, z, b) {
    (t(x) %*% sinv) %*% (y - x %*% beta - z %*% b)
  },
  x = Xi, sinv = Sigmai.inv, y = yi, z = Zi, b = Eb,
  SIMPLIFY = TRUE)

  rownames(sbeta) <- names(beta)

  #-----------------------------------------------------

  # D

  Dinv <- solve(D)

  D.inds <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  dimnames(D.inds) <- NULL

  delta.D <- lapply(1:nrow(D.inds), function(x, ind) {
    mat <- matrix(0, nrow = nrow(D), ncol = ncol(D))
    ii <- ind[x, , drop = FALSE]
    mat[ii[1], ii[2]] <- mat[ii[2], ii[1]] <- 1
    mat
  }, ind = D.inds[, 2:1, drop = FALSE])

  term1 <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })

  sDi <- function(i) {
    mapply(function(b, pb) {
      out <- 0.5 * crossprod(b, b * pb) %*% (Dinv %*% delta.D[[i]] %*% Dinv) / nrow(b)
      term1[i] + sum(diag(out))
    },
    b = bi.y, pb = pb.yt,
    SIMPLIFY = TRUE)
  }

  sD <- sapply(1:nrow(D.inds), sDi)
  sD <- t(sD)
  rownames(sD) <- paste0("D", D.inds[, 1], ",", D.inds[, 2])

  #-----------------------------------------------------

  # gamma

  sgamma <- gammaUpdate_approx(bi.y, Zit.fail, expvstargam, pb.yt, haz.hat,
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
      resids <- t(residFixed) %*% (residFixed - 2*(z %*% b.k)) + sum(diag(crossprod(z) %*% bbT.k))
      (-0.5 * nik[k] / sigma2[k]) + (0.5 * resids / sigma2[k]^2)
    },
    y = yik[[k]], x = Xik.list[[k]], z = Zik.list[[k]], b = Eb, b2 = EbbT, nik = nik)
  }

  rownames(ssigma2) <- paste0("sigma2_", 1:K)

  #-----------------------------------------------------

  si <- rbind(sD, sbeta, ssigma2, sgamma)

  H <- matrix(0, nrow(si), nrow(si))
  for (j in 1:ncol(si)) {
    H <- H + tcrossprod(si[, j])
  }
  # Although RHS term = 0 in theory, in practice with MC integration
  # not all terms are vanishingly small, so we add it in
  H <- H - (rowSums(si) %*% t(rowSums(si))) / ncol(si)
  rownames(H) <- colnames(H) <- rownames(si)

  return(H)

}
