# Hazard contributions for event times with indices after 'start' and up to
# data$tj.ind. Subsetting first avoids subtracting large cumulative hazards.
.mjoint_survival_terms <- function(theta, data, start = 0L) {
  r <- sum(data$r)
  j <- data$tj.ind
  if (start > j) stop("Prediction time precedes the conditioning time.")
  if (start == j) {
    return(list(A = matrix(numeric(), 0L, r), logh = numeric()))
  }
  indices <- seq.int(start + 1L, j)
  haz <- theta$haz[indices]
  if (any(!is.finite(haz)) || any(haz < 0)) {
    stop("Baseline-hazard increments must be finite and nonnegative.")
  }
  if (length(theta$gamma) != data$q + data$K || any(!is.finite(theta$gamma))) {
    stop("Survival coefficients must be finite and match the fitted model.")
  }
  # Zero increments contribute nothing, even for extreme linear predictors.
  indices <- indices[haz > 0]
  logh <- log(haz[haz > 0])
  A <- (data$IW.fail %*% data$Z.fail)[indices, , drop = FALSE]
  A <- sweep(A, 2, rep(utils::tail(theta$gamma, data$K), data$r), "*")
  if (data$q > 0) logh <- logh + sum(data$v * utils::head(theta$gamma, data$q))
  if (any(!is.finite(c(A, logh)))) {
    stop("Survival contributions are not finite.")
  }
  list(A = A, logh = logh)
}

#' @keywords internal
S <- function(b, theta, data, start = 0L) {
  if (length(b) != sum(data$r) || any(!is.finite(b))) {
    stop("Random effects must be finite and match the fitted model.")
  }
  terms <- .mjoint_survival_terms(theta, data, start)
  log.haz <- as.vector(terms$logh + terms$A %*% b)
  if (anyNA(log.haz)) stop("The conditional survival probability is undefined.")
  exp(-sum(exp(log.haz)))
}
