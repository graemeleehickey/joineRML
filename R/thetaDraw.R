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

  # Sample residual variances on the log scale, retaining their joint
  # covariance with the other parameters through the delta method.
  sigma2 <- theta.mean[["sigma2"]]
  if (length(sigma2) != object$dims$K ||
      any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    stop("Simulated predictions require finite, strictly positive fitted residual variances.")
  }
  draw.mean <- unlist(as.relistable(theta.mean))
  sigma.inds <- which(rep(names(theta.mean), lengths(theta.mean)) == "sigma2")
  if (any(!is.finite(draw.mean))) {
    stop("Simulated predictions require finite fitted parameters.")
  }
  if (!is.matrix(theta.var) || !is.numeric(theta.var) ||
      !identical(dim(theta.var), rep(length(draw.mean), 2L)) ||
      any(!is.finite(theta.var)) ||
      !isSymmetric(theta.var, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("Simulated predictions require a finite, symmetric covariance matrix matching the fitted parameters.")
  }
  jacobian <- rep(1, length(draw.mean))
  jacobian[sigma.inds] <- 1 / sigma2
  draw.mean[sigma.inds] <- log(sigma2)
  draw.var <- sweep(sweep(theta.var, 1, jacobian, "*"), 2, jacobian, "*")
  if (any(!is.finite(draw.var))) {
    stop("The covariance matrix for simulated predictions is not finite after transforming residual variances to the log scale.")
  }
  theta.samp <- withCallingHandlers(
    mvtnorm::rmvnorm(n = 1, mean = draw.mean, sigma = draw.var),
    warning = function(w) {
      stop(paste("Cannot simulate parameters from the fitted covariance matrix:",
                 conditionMessage(w)), call. = FALSE)
    }
  )
  theta.samp[, sigma.inds] <- exp(theta.samp[, sigma.inds])
  if (any(!is.finite(theta.samp)) || any(theta.samp[, sigma.inds] <= 0)) {
    stop("Simulated parameters are not finite or residual variances are not strictly positive after transformation; check the fitted variances and their uncertainty.")
  }
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
