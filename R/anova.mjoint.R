#' Compare fitted joint models
#'
#' @description Compare two or more fitted joint models using their joint
#'   log-likelihoods, AIC, and BIC. For nested covariate models, also calculate
#'   likelihood-ratio statistics and approximate chi-squared p-values.
#'
#' @param object an object inheriting from class \code{mjoint}.
#' @param ... additional fitted \code{mjoint} objects, ordered from the reduced
#'   model to the full model when \code{test = TRUE}.
#' @param test logical: calculate likelihood-ratio tests between successive
#'   models? Default is \code{TRUE}. With \code{FALSE}, report the stored
#'   log-likelihoods, AIC, and BIC without testing or further integration.
#' @param nMC an integer of at least 2 giving the number of Halton integration
#'   points per subject used to re-evaluate each likelihood when testing.
#'   Default is \code{10000}. The first \code{floor(nMC / 2)} points are also
#'   used to assess numerical stability.
#' @param tol a positive finite tolerance for the absolute change in the
#'   likelihood-ratio statistic between the two integration sizes. Default is
#'   \code{0.1}. This is a numerical diagnostic, not a significance threshold
#'   or a bound on integration error.
#'
#' @details All models must have converged and use the same subjects,
#'   longitudinal responses and measurement times, and survival responses.
#'   Longitudinal response expressions must also agree, so comparisons between
#'   different response transformations are rejected. Models fitted with
#'   \code{pfs = FALSE} are supported. Older fitted objects that lack the stored
#'   fitting matrices must be refitted with the current version of \code{mjoint}.
#'
#'   Tests require unchanged random-effects and association structures, and
#'   nested longitudinal fixed-effects and survival-covariate design matrices.
#'   Each successive model must add at least one identifiable coefficient.
#'   Nesting is checked using matrix column spaces, allowing different factor
#'   contrasts. The reference distribution is chi-squared, with degrees of
#'   freedom equal to the difference in parameter counts. The common baseline
#'   hazard has the same event-time support in both models and does not add to
#'   this difference. Tests that change the random-effects structure, including
#'   tests of zero variance components, are not supported. Non-nested models on
#'   the same observations can be compared with \code{test = FALSE}.
#'
#'   When testing, the observed joint likelihood is re-evaluated at the fitted
#'   coefficients using common deterministic Halton points for the conditional
#'   Gaussian random effects. No models are refitted. The table's log-likelihood,
#'   AIC, and BIC use these re-evaluated values, so they may differ from the
#'   corresponding values stored in the fitted objects. BIC uses the number of
#'   subjects. Neither the fitted objects nor R's random-number state is changed.
#'
#'   A warning is issued and the p-value is set to \code{NA} if a likelihood-ratio
#'   statistic is negative or changes by more than \code{tol} when the number of
#'   integration points increases. Increase \code{nMC} to investigate integration
#'   error. A negative statistic can also indicate that a fitted model needs
#'   tighter convergence settings or better initial values. Passing this
#'   diagnostic does not establish that either model reached the maximum, or
#'   guarantee accurate integration or validity of the asymptotic approximation.
#'
#'   A single-model call is not supported; no sequential term deletion or
#'   automatic model fitting is performed.
#'
#' @return An object of class \code{anova}, also inheriting from
#'   \code{data.frame}, with one row per model. Columns \code{df}, \code{logLik},
#'   \code{AIC}, and \code{BIC} describe each model. When testing, \code{L.Ratio},
#'   \code{Df}, and \code{Pr(>Chisq)} describe comparisons with the preceding
#'   model, and \code{LR.change} is the absolute change in the statistic between
#'   the two integration sizes. The \code{integration} attribute records the
#'   method, integration sizes, tolerance, and likelihoods at the smaller size.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @seealso \code{\link{mjoint}}, \code{\link{logLik.mjoint}}.
#' @export
#'
#' @examples
#' \dontrun{
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.lvmi), ]
#' reduced <- mjoint(
#'   formLongFixed = log.lvmi ~ time,
#'   formLongRandom = ~ 1 | num,
#'   formSurv = Surv(fuyrs, status) ~ age,
#'   data = hvd, timeVar = "time", seed = 123)
#' full <- update(reduced, formLongFixed = log.lvmi ~ time + sex)
#' anova(reduced, full)
#' # Check sensitivity to more integration points
#' anova(reduced, full, nMC = 20000)
#' # Descriptive comparison using the stored likelihoods
#' anova(reduced, full, test = FALSE)
#' }
anova.mjoint <- function(object, ..., test = TRUE, nMC = 10000L, tol = 0.1) {
  models <- list(object, ...)
  if (length(models) < 2L) {
    stop("Supply at least two fitted mjoint models to compare.")
  }
  if (!all(vapply(models, inherits, logical(1), what = "mjoint"))) {
    stop("All objects must be fitted mjoint models.")
  }
  if (!is.logical(test) || length(test) != 1L || is.na(test)) {
    stop("'test' must be TRUE or FALSE.")
  }
  if (!is.numeric(nMC) || is.complex(nMC) || length(nMC) != 1L ||
      !is.finite(nMC) || nMC < 2 || nMC > .Machine$integer.max ||
      nMC != trunc(nMC)) {
    stop("'nMC' must be a single integer between 2 and .Machine$integer.max.")
  }
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive finite number.")
  }
  if (!all(vapply(models, function(x) isTRUE(x$conv), logical(1)))) {
    stop("All models must have converged before comparison.")
  }
  if (any(vapply(models, function(x) is.null(x$dmats), logical(1)))) {
    stop("A model lacks stored fitting matrices; refit it with the current version of mjoint().")
  }

  labels <- vapply(as.list(substitute(list(object, ...)))[-1L], deparse1,
                   character(1))
  if (!is.null(names(models))) {
    named <- nzchar(names(models))
    labels[named] <- names(models)[named]
  }
  labels <- make.unique(labels)
  observed <- lapply(models, mjoint_observed_data)
  for (i in seq.int(2L, length(models))) {
    if (!identical(observed[[1]], observed[[i]])) {
      stop("Models must use the same subjects, responses, response transformations, and measurement times.")
    }
  }

  ll <- vapply(models, function(x) as.numeric(logLik(x)), numeric(1))
  df <- vapply(models, function(x) attr(logLik(x), "df"), numeric(1))
  n <- vapply(models, function(x) attr(logLik(x), "nobs"), numeric(1))
  if (any(!is.finite(c(ll, df, n))) || any(n <= 0)) {
    stop("Models must have finite log-likelihoods and valid parameter and subject counts.")
  }

  if (test) {
    for (i in seq.int(2L, length(models))) {
      reduced <- models[[i - 1L]]
      full <- models[[i]]
      if (!mjoint_same_random_structure(reduced, full)) {
        stop("Likelihood-ratio tests require the same random-effects and association structures; use test = FALSE for a descriptive comparison.")
      }
      x0 <- do.call(rbind, reduced$dmats$l$Xi)
      x1 <- do.call(rbind, full$dmats$l$Xi)
      # A constant survival predictor is absorbed by the baseline hazard.
      v0 <- cbind(1, reduced$sfit$x)
      v1 <- cbind(1, full$sfit$x)
      if (!mjoint_nested_matrix(x0, x1) || !mjoint_nested_matrix(v0, v1)) {
        stop("Covariate models must be identifiable and nested in the supplied order; use test = FALSE for non-nested comparisons.")
      }
    }
    if (any(diff(df) <= 0)) {
      stop("Each successive model must add at least one identifiable coefficient.")
    }
    nodes <- matrix(randtoolbox::halton(nMC, dim = sum(object$dims$r),
                                       normal = TRUE), nrow = nMC)
    likelihoods <- vapply(models, mjoint_comparison_loglik, numeric(2), nodes = nodes)
    if (any(!is.finite(likelihoods))) {
      stop("Likelihood re-evaluation was non-finite; check the fitted models and integration settings.")
    }
    ll <- likelihoods[2L, ]
    lr <- 2 * diff(ll)
    lr.change <- abs(lr - 2 * diff(likelihoods[1L, ]))
    reliable <- lr >= 0 & lr.change <= tol
    pvalue <- rep(NA_real_, length(lr))
    pvalue[reliable] <- pchisq(lr[reliable], diff(df)[reliable], lower.tail = FALSE)
    if (any(lr < 0)) {
      warning("A likelihood-ratio statistic is negative; its p-value is withheld. Check integration precision and refit with tighter convergence settings or better initial values.", call. = FALSE)
    }
    if (any(lr.change > tol)) {
      warning("A likelihood-ratio statistic is unstable between integration sizes; its p-value is withheld. Increase nMC and check LR.change.", call. = FALSE)
    }
  }

  out <- data.frame(df = df, logLik = ll, AIC = -2 * ll + 2 * df,
                    BIC = -2 * ll + log(n) * df, row.names = labels)
  heading <- "Joint model comparison (maximum likelihood)\n"
  if (test) {
    out$L.Ratio <- c(NA_real_, lr)
    out$Df <- c(NA_real_, diff(df))
    out[["Pr(>Chisq)"]] <- c(NA_real_, pvalue)
    out$LR.change <- c(NA_real_, lr.change)
    attr(out, "integration") <- list(method = "halton", nMC = nMC,
      nMC.small = floor(nMC / 2), tol = tol, logLik.small = likelihoods[1L, ])
    heading <- c(heading, paste0("Likelihoods re-evaluated with ", nMC,
      " common Halton points per subject. LR.change compares ",
      floor(nMC / 2), " and ", nMC, " points (tolerance ", tol, ").\n"))
  } else {
    heading <- c(heading, "Stored likelihoods; no likelihood-ratio tests.\n")
  }
  attr(out, "heading") <- heading
  class(out) <- c("anova", "data.frame")
  out
}

# Compare the observations actually used, independently of unused data columns.
mjoint_observed_data <- function(object) {
  longitudinal <- lapply(seq_len(object$dims$K), function(k) {
    d <- object$data[[k]]
    list(response = deparse1(object$formLongFixed[[k]][[2L]]),
         id = as.character(d[[object$id]]),
         time = as.numeric(d[[object$timeVar[k]]]))
  })
  surv <- object$dmats$t$survdat2
  list(longitudinal = longitudinal,
       y = lapply(object$dmats$l$yi, function(y) as.numeric(y)),
       id = as.character(surv$id), time = as.numeric(surv$T),
       status = as.numeric(surv$delta))
}

mjoint_same_random_structure <- function(x, y) {
  same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 0,
                                         check.attributes = FALSE))
  same(x$dims$r, y$dims$r) &&
    same(x$dmats$l$Zi, y$dmats$l$Zi) &&
    same(x$dmats$z$Zi.fail, y$dmats$z$Zi.fail) &&
    same(x$dmats$z$IW.fail, y$dmats$z$IW.fail)
}

mjoint_nested_matrix <- function(reduced, full) {
  if (nrow(reduced) != nrow(full) || any(!is.finite(c(reduced, full)))) {
    return(FALSE)
  }
  rank <- function(x) {
    norms <- sqrt(colSums(x^2))
    x <- sweep(x, 2L, ifelse(norms > 0, norms, 1), "/")
    qr(x, tol = 1e-7)$rank
  }
  rank(reduced) == ncol(reduced) && rank(full) == ncol(full) &&
    rank(cbind(full, reduced)) == ncol(full)
}

# Evaluate one subject at a time to avoid storing samples for the whole cohort.
# Log-scale averaging also avoids underflow in small survival likelihoods.
mjoint_comparison_loglik <- function(object, nodes) {
  theta <- object$coefficients
  l <- object$dmats$l
  t <- object$dmats$t
  z <- object$dmats$z
  Dinv <- solve(theta$D)
  association <- rep(utils::tail(theta$gamma, l$K), l$r)
  half <- seq_len(floor(nrow(nodes) / 2))
  totals <- c(0, 0)
  logmean <- function(x) {
    largest <- max(x)
    if (is.infinite(largest)) return(largest)
    largest + log(mean(exp(x - largest)))
  }
  for (i in seq_len(l$n)) {
    Z <- l$Zi[[i]]
    mu <- as.vector(l$Xi[[i]] %*% theta$beta)
    variance <- rep(theta$sigma2, l$nik[[i]])
    A <- solve(Dinv + crossprod(Z, Z / variance))
    m <- as.vector(A %*% crossprod(Z, (l$yi[[i]] - mu) / variance))
    b <- sweep(nodes %*% chol(A), 2L, m, "+")
    fy <- mvtnorm::dmvnorm(as.vector(l$yi[[i]]), mu,
      sigma = Z %*% theta$D %*% t(Z) + diag(variance, nrow = length(variance)),
      log = TRUE)
    h <- t$survdat2.list[[i]]
    nt <- h$tj.ind
    log.surv <- rep(0, nrow(nodes))
    if (nt > 0) {
      design <- t(z$IW.fail[[i]] %*% z$Zi.fail[[i]])
      eta <- b %*% (design * association)
      if (t$q > 0) {
        eta <- eta + sum(t$V[[i]] * theta$gamma[seq_len(t$q)])
      }
      log.haz <- log(theta$haz[seq_len(nt)])
      log.surv <- -rowSums(exp(sweep(eta, 2L, log.haz, "+")))
      if (h$delta == 1) {
        log.surv <- log.surv + log.haz[nt] + eta[, nt]
      }
    }
    totals <- totals + fy + c(logmean(log.surv[half]), logmean(log.surv))
  }
  totals
}
