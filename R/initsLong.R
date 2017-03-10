#' @keywords internal
initsLong <- function(lfit, inits, l, z, K, p, r, tol.em, verbose) {

  D <- Matrix::bdiag(lapply(lfit, function(u) matrix(nlme::getVarCov(u),
                                                     dim(nlme::getVarCov(u)))))
  D <- as.matrix(D)
  D.names <- c()
  for (k in 1:K) {
    D.names.k <- paste0(rownames(nlme::getVarCov(lfit[[k]])), "_", k)
    D.names <- c(D.names, D.names.k)
  }
  rownames(D) <- colnames(D) <- D.names

  beta <- do.call("c", lapply(lfit, fixef))
  names(beta) <- paste0(names(beta), "_", rep(1:K, p))

  sigma2 <- unlist(lapply(lfit, function(u) u$sigma))^2

  if ((K > 1) && !all(c("beta", "D", "sigma2") %in% names(inits))) {
    message("Running multivariate LMM EM algorithm to establish initial parameters...")
    out <- mvlme(thetaLong = list("beta" = beta, "D" = D, "sigma2" = sigma2),
                 l = l, z = z, tol.em = tol.em, verbose = verbose)
    message("Finished multivariate LMM EM algorithm...")
  } else {
    out <- list("D" = D, "beta" = beta, "sigma2" = sigma2)
  }

  # over-ride with user-specified inits
  if ("beta" %in% names(inits)) {
    if (length(inits$beta) != sum(p)) {
      stop("Dimension of beta inits does not match model.")
    }
    beta <- inits$beta
    names(beta) <- names(out[["beta"]])
    out[["beta"]] <- beta
  }
  if ("D" %in% names(inits)) {
    if (nrow(inits$D) != sum(r)) {
      stop("Dimension of D inits does not match model.")
    }
    is.posdef <- all(eigen(inits$D)$values > 0)
    if (is.posdef) {
      D <- inits$D
      rownames(D) <- colnames(D) <- rownames(out[["D"]])
      out[["D"]] <- D
    } else {
      warning("Initial parameter matrix D is non positive definite: falling back to automated value")
    }
  }
  if ("sigma2" %in% names(inits)) {
    sigma2 <- inits$sigma2
    if (length(sigma2) != K) {
      stop("Dimension of sigma2 inits does not match model.")
    }
    names(sigma2) <- paste0("sigma2_", 1:K)
    out[["sigma2"]] <- sigma2
  }

  return(out)

}
