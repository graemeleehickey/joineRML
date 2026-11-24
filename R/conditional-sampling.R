# Gaussian longitudinal conditioning followed by exact rejection sampling of
# the survival contribution. All calculations hold the parameter draw fixed.
.mjoint_conditional_base <- function(theta, data) {
  r <- sum(data$r)
  D <- as.matrix(theta$D)
  if (!all(dim(D) == c(r, r)) || any(!is.finite(D)) ||
      !isSymmetric(D, check.attributes = FALSE)) {
    stop("The random-effects covariance must be finite and symmetric.")
  }
  if (length(theta$sigma2) != data$K || any(!is.finite(theta$sigma2)) ||
      any(theta$sigma2 <= 0)) {
    stop("Conditional sampling requires finite, strictly positive residual variances.")
  }
  LD <- tryCatch(t(chol(D)), error = function(e) {
    stop("The random-effects covariance must be positive definite.")
  })
  residual.sd <- sqrt(rep(theta$sigma2, data$nk))
  W <- (data$Z %*% LD) / residual.sd
  residual <- as.vector(data$y - data$X %*% theta$beta) / residual.sd
  precision <- diag(r) + crossprod(W)
  if (any(!is.finite(c(precision, residual)))) {
    stop("The longitudinal conditional distribution is not finite.")
  }
  C <- chol2inv(chol(precision))
  m <- as.vector(LD %*% C %*% crossprod(W, residual))
  L <- LD %*% t(chol(C))
  terms <- .mjoint_survival_terms(theta, data)
  # Survival terms that do not depend on b cancel from its conditional density.
  keep <- rowSums(abs(terms$A)) > 0
  A <- terms$A[keep, , drop = FALSE]
  B <- A %*% L
  logh <- as.vector(terms$logh[keep] + A %*% m)
  if (any(!is.finite(c(m, L, B, logh)))) {
    stop("The conditional random-effects distribution is not finite.")
  }
  list(m = m, L = L, B = B, logh = logh)
}

.mjoint_rejection_envelope <- function(base) {
  B <- base$B
  logh <- base$logh
  r <- length(base$m)
  objective <- function(z) sum(z^2) / 2 + sum(exp(logh + B %*% z))
  gradient <- function(z) as.vector(z + crossprod(B, exp(logh + B %*% z)))
  hessian <- function(z) diag(r) + crossprod(B, B * as.vector(exp(logh + B %*% z)))
  anchor <- rep(0, r)
  if (nrow(B)) {
    if (!is.finite(objective(anchor))) {
      stop("Cannot construct a finite conditional-sampling envelope; check the fitted model and parameter uncertainty.")
    }
    opt <- stats::nlminb(anchor, objective, gradient, hessian,
                        control = list(iter.max = 100, eval.max = 300,
                                       rel.tol = 1e-11, x.tol = 1e-9))
    anchor <- opt$par
  }
  logw <- as.vector(logh + B %*% anchor)
  # The tangent bound is valid at any finite anchor, including an approximate
  # mode. Its proposal mean is the negative survival gradient, not the anchor.
  center <- -as.vector(crossprod(B, exp(logw)))
  if (any(!is.finite(c(anchor, logw, center)))) {
    stop("The conditional-sampling envelope is not finite.")
  }
  c(base, list(anchor = anchor, center = center, logw = logw))
}

# log(exp(x) - 1 - x), evaluated without cancellation or exponential overflow.
.mjoint_log_tangent_gap <- function(x) {
  if (anyNA(x)) stop("The conditional-sampling acceptance probability is undefined.")
  out <- rep(Inf, length(x))
  small <- is.finite(x) & abs(x) < 1e-4
  out[small] <- 2 * log(abs(x[small])) - log(2) +
    log1p(x[small] / 3 + x[small]^2 / 12 + x[small]^3 / 60 + x[small]^4 / 360)
  large <- is.finite(x) & x > 50
  out[large] <- x[large] + log1p(-(1 + x[large]) * exp(-x[large]))
  middle <- is.finite(x) & !small & !large
  out[middle] <- log(expm1(x[middle]) - x[middle])
  out
}

.mjoint_rejection_log_accept <- function(z, envelope) {
  delta <- as.vector(envelope$B %*% (z - envelope$anchor))
  value <- -sum(exp(envelope$logw + .mjoint_log_tangent_gap(delta)))
  if (is.na(value) || value > 0) {
    stop("The conditional-sampling acceptance probability is invalid.")
  }
  value
}

.mjoint_rejection_draw <- function(envelope, max.proposals) {
  for (proposal in seq_len(max.proposals)) {
    z <- envelope$center + rnorm(length(envelope$m))
    if (log(runif(1)) <= .mjoint_rejection_log_accept(z, envelope)) {
      b <- as.vector(envelope$m + envelope$L %*% z)
      if (any(!is.finite(b))) stop("The sampled random effects are not finite.")
      return(list(b = b, proposals = proposal))
    }
  }
  stop(paste0("Conditional sampling exhausted max.proposals = ", max.proposals,
              "; no parameter draw was discarded. Check the fit and its uncertainty, ",
              "or increase max.proposals if the conditional distribution is difficult to sample."))
}

.mjoint_prediction_draw <- function(object, data, max.proposals, iteration, M) {
  tryCatch({
    theta <- thetaDraw(object)
    envelope <- .mjoint_rejection_envelope(.mjoint_conditional_base(theta, data))
    draw <- .mjoint_rejection_draw(envelope, max.proposals)
    c(list(theta = theta), draw)
  }, error = function(e) {
    stop(sprintf("Simulation %s of %s at landmark %s failed: %s",
                 iteration, M, format(data$tobs), conditionMessage(e)), call. = FALSE)
  })
}

.mjoint_simulation_controls <- function(M, max.proposals, ci) {
  for (name in c("M", "max.proposals")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(paste0("'", name, "' must be a finite, positive integer."))
    }
  }
  if (!is.numeric(ci) || length(ci) != 1L || !is.finite(ci) || ci < 0 || ci > 1) {
    stop("ci argument has been misspecified.")
  }
}
