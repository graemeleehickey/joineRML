#' Internal function for generating initial parameters for the longitudinal
#' sub-model
#'
#' @keywords internal
initsLong <- function(lfit, inits, l, z, K, p, tol.em, verbose) {

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
    beta <- inits$beta
    names(beta) <- names(out[["beta"]])
    out[["beta"]] <- beta
  }
  if ("D" %in% names(inits)) {
    D <- inits$D
    rownames(D) <- colnames(D) <- rownames(out[["D"]])
    out[["D"]] <- D
  }
  if ("sigma2" %in% names(inits)) {
    sigma2 <- inits$sigma2
    names(sigma2) <- paste0("sigma2_", 1:K)
    out[["sigma2"]] <- sigma2
  }

  return(out)

}
