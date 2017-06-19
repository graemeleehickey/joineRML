#' @keywords internal
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils relist as.relistable
thetaDraw <- function(object) {

  # Mean
  theta.mean <- coef(object)
  theta.mean <- theta.mean[-which(names(theta.mean) == "haz")]
  D.inds <- which(lower.tri(theta.mean[["D"]], diag = TRUE), arr.ind = TRUE)
  theta.mean[["D"]] <- theta.mean[["D"]][D.inds]

  # Variance
  theta.var <- vcov(object)

  theta.samp <- mvtnorm::rmvnorm(n = 1,
                                 mean = unlist(as.relistable(theta.mean)),
                                 sigma = theta.var)
  theta.samp <- relist(theta.samp, skeleton = theta.mean)
  D <- matrix(0, nrow = max(D.inds), ncol = max(D.inds))
  D[D.inds] <- theta.samp[["D"]]
  if (sum(object$dims$r) > 1) {
    D <- D + t(D) - diag(diag(D))
    D <- Matrix::nearPD(D)
    theta.samp[["D"]] <- as.matrix(D$mat)
  } else {
    D <- max(D, 1e-08)
    theta.samp[["D"]] <- D
  }

  # Baseline hazard
  haz <- baseHaz(object, se = TRUE)
  haz.samp <- rnorm(nrow(haz), mean = haz$haz, sd = haz$se)
  haz.samp <- pmax(haz.samp, min(1e-06, haz$haz))

  theta.samp[["haz"]] <- haz.samp
  return(theta.samp)

}
