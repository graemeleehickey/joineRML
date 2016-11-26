#' Standard errors via bootstrap for an \code{mjoint} object
#'
#' This function takes a model fit from an \code{mjoint} object and calculates
#' standard errors and confidence intervals for the main longitudinal and
#' survival coefficient parameters, including the latent association parameters,
#' using bootstrapping (Efron and Tibshirani, 2000).
#'
#' @param object an object inheriting from class \code{mjoint} for a joint model
#'   of time-to-event and multivariate longitudinal data.
#' @param nboot the number of bootstrap samples. Default is \code{nboot=100}.
#' @param ci the confidence interval to be estimated using the
#'   percentile-method. Default is \code{ci=0.95} for a 95\% confidence
#'   interval.
#' @param use.mle logical: should the algorithm use the maximiser from the
#'   converged model in \code{object} as initial values for coefficients in each
#'   bootstrap iteration. Default is \code{use.mle=TRUE}.
#' @param progress logical: should a progress bar be shown on the console to
#'   indicate the percentage of bootstrap iterations completed? Default is
#'   \code{progress=TRUE}.
#' @inheritParams mjoint
#'
#' @details Standard errors and confidence intervals are obtained by repeated
#'   fitting of the requisite joint model to bootstrap samples of the original
#'   longitudinal and time-to-event data. Note that bootstrap is done by
#'   sampling subjects, not individual records.
#'
#' @note This function is computationally intensive. A dynamic progress bar is
#'   displayed showing the percentage of bootstrap models fitted.
#'
#'   Due to random sampling, an \code{mjoint} model fitted to some bootstrap
#'   samples may not converge within the specified control parameter settings.
#'   The \code{bootSE} code discards any models that failed to converge when
#'   calculating the standard errors and confidence intervals. If a large
#'   proportion of models have failed to converge, it is likely that it will
#'   need to be refitted with changes to the \code{control} arguments.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords multivariate survival methods
#' @seealso \code{\link{mjoint}} for approximate standard errors.
#'
#' @references
#'
#' Efron B, Tibshirani R. \emph{An Introduction to the Bootstrap.} 2000; Boca
#' Raton, FL: Chapman & Hall/CRC.
#'
#' @return An object of class \code{bootSE}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a joint model with bivariate longitudinal outcomes
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#'
#' fit2 <- mjoint(
#'     formLongFixed = list("grad" = log.grad ~ time + sex + hs,
#'                          "lvmi" = log.lvmi ~ time + sex),
#'     formLongRandom = list("grad" = ~ 1 | num,
#'                           "lvmi" = ~ time | num),
#'     formSurv = Surv(fuyrs, status) ~ age,
#'     data = list(hvd, hvd),
#'     inits = list("gamma" = c(0.11, 1.51, 0.80)),
#'     timeVar = "time",
#'     verbose = TRUE)
#'
#' fit2.boot <- bootSE(fit1, 50, use.mle = TRUE, control = list(
#'     earlyPhase = 25, convCrit = "either",
#'     tol0 = 6e-03, tol2 = 6e-03, mcmaxIter = 60))
#' }
bootSE <- function(object, nboot = 100, ci = 0.95, use.mle = TRUE,
                   verbose = FALSE, control = list(), progress = TRUE,
                   ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  t.start <- Sys.time()

  # Extract from fitted model
  formLongFixed <- object$formLongFixed
  formLongRandom <- object$formLongRandom
  formSurv <- object$formSurv
  timeVar <- object$timeVar
  K <- object$dims$K

  # Use fitted model MLE as initial values?
  if (use.mle) {
    theta.hat <- object$coefficients
    theta.hat <- theta.hat[-which(names(theta.hat) == "haz")]
  } else {
    theta.hat <- NULL
  }

  out <- list()
  conv.status <- vector(length = nboot)
  if (progress) {
    cat("\n\n")
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
  }
  for (b in 1:nboot) {
    # bootstrap sample data
    data.boot <- sampleData(object = object)
    # fit joint model
    suppressMessages(
      fit.boot <- mjoint(formLongFixed = formLongFixed,
                         formLongRandom = formLongRandom,
                         formSurv = formSurv,
                         data = data.boot$longData.boot,
                         survData = data.boot$survData.boot,
                         timeVar = timeVar,
                         inits = theta.hat,
                         verbose = verbose,
                         se.approx = FALSE,
                         postRE = FALSE,
                         control = control,
                         ...)
    )
    out[[b]] <- fit.boot$coefficients
    conv.status[b] <- fit.boot$conv
    if (progress) {
      utils::setTxtProgressBar(pb, b)
    }
  }
  if (progress) {
    close(pb)
  }

  # Checks
  if (mean(conv.status) <= 0.1) {
    stop("Cannot estimate SEs: less than 10% of bootstrap models converged.")
  }
  out <- out[conv.status]

  # Bootstrap estimates
  beta.mat   <- sapply(out, function(u) u$beta)
  gamma.mat  <- sapply(out, function(u) u$gamma)
  sigma2.mat <- sapply(out, function(u) u$sigma2)

  # If K = 1, sigma2 might be a vector, so need to coerce to matrix
  vec2mat <- function(x) {
    if (!is.matrix(x) && is.vector(x)) {
      cnam <- names(x)[1]
      x <- matrix(x, nrow = 1)
      rownames(x) <- cnam
    }
    return(x)
  }
  sigma2.mat <- vec2mat(sigma2.mat)

  # CI function
  qFun <- function(x) {
    pr <- (1 - ci) / 2
    quantile(x, c(pr, 1 - pr))
  }

  beta.ci <- apply(beta.mat, 1, qFun)
  gamma.ci <- apply(gamma.mat, 1, qFun)
  sigma2.ci <- apply(sigma2.mat, 1, qFun)

  sigma2.se <- apply(sigma2.mat, 1, sd)
  beta.se <- apply(beta.mat, 1, sd)
  gamma.se <- apply(gamma.mat, 1, sd)

  t.total <- Sys.time() - t.start

  out.list <- list(beta.ci = beta.ci, gamma.ci = gamma.ci, sigma2.ci = sigma2.ci)
  out.list$beta.se <- beta.se
  out.list$gamma.se <- gamma.se
  out.list$sigma2.se <- sigma2.se
  out.list$boot.time <- t.total
  out.list$nboot <- nboot
  out.list$ci <- ci
  out.list$coefficients <- object$coefficients
  out.list$conv <- sum(conv.status)

  class(out.list) <- "bootSE"
  invisible(out.list)

}
