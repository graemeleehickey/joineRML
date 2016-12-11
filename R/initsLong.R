#' Internal function for generating initial parameters for the longitudinal
#' sub-model
#'
#' @keywords internal
initsLong <- function(lfit, K, p) {

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

  return(list("D" = D, "beta" = beta, "sigma2" = sigma2))

}
