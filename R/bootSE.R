#' Standard errors via bootstrap for an \code{mjoint} object
#'
#' @description This function takes a model fit from an \code{mjoint} object and
#'   calculates standard errors and confidence intervals for the main
#'   longitudinal and survival coefficient parameters, including the latent
#'   association parameters, using bootstrapping (Efron and Tibshirani, 2000).
#'
#' @param object an object inheriting from class \code{mjoint} for a joint model
#'   of time-to-event and multivariate longitudinal data.
#' @param nboot the number of bootstrap samples. Default is \code{nboot=100}.
#' @param ci the confidence interval to be estimated using the
#'   percentile-method. Default is \code{ci=0.95} for a 95\% confidence
#'   interval.
#' @param use.mle logical: should the algorithm use the maximizer from the
#'   converged model in \code{object} as initial values for coefficients in each
#'   bootstrap iteration. Default is \code{use.mle=TRUE}.
#' @param progress logical: should a progress bar be shown on the console to
#'   indicate the percentage of bootstrap iterations completed? Default is
#'   \code{progress=TRUE}.
#' @param ncores integer: if more than one core is available, then the \code{bootSE}
#'   function can run in parallel via the \code{\link[foreach]{foreach}}
#'   function. By default, \code{ncores=1}, which defaults to serial mode. Note
#'   that if \code{ncores}>1, then \code{progress} is set to \code{FALSE} by
#'   default, as it is not possible to display progress bars for parallel
#'   processes at the current time.
#' @param safe.boot logical: should each bootstrap replication be wrapped in a
#'   \code{\link[base]{tryCatch}} statement to catch errors (e.g. during the
#'   optimisation progress)? When model fitting throws errors, a new bootstrap
#'   sample is drawn for the current iteration and the model is re-fit; this
#'   process continuex until a model fits successfully. Default is \code{FALSE}.
#' @inheritParams mjoint
#'
#' @details Standard errors and confidence intervals are obtained by repeated
#'   fitting of the requisite joint model to bootstrap samples of the original
#'   longitudinal and time-to-event data. Note that bootstrap is done by
#'   sampling subjects, not individual records.
#'
#' @note This function is computationally intensive. A dynamic progress bar is
#'   displayed showing the percentage of bootstrap models fitted. On computer
#'   systems with more than one core available, computational time can be
#'   reduced by passing the argument \code{ncores} (with integer value >1) to
#'   \code{bootSE}, which implements parallel processing via the
#'   \code{\link[foreach]{foreach}} package. \strong{Note:} if parallel
#'   processing is implemented, then the progress bar is not displayed.
#'
#'   Due to random sampling, an \code{mjoint} model fitted to some bootstrap
#'   samples may not converge within the specified control parameter settings.
#'   The \code{bootSE} code discards any models that failed to converge when
#'   calculating the standard errors and confidence intervals. If a large
#'   proportion of models have failed to converge, it is likely that it will
#'   need to be refitted with changes to the \code{control} arguments.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords multivariate survival methods
#' @seealso \code{\link{mjoint}} for approximate standard errors.
#'
#' @references
#'
#' Efron B, Tibshirani R. \emph{An Introduction to the Bootstrap.} 2000; Boca
#' Raton, FL: Chapman & Hall/CRC.
#'
#' @return An object of class \code{bootSE}.
#' @import foreach
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a joint model with bivariate longitudinal outcomes
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#'
#' fit <- mjoint(
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
#' fit.boot <- bootSE(fit, 50, use.mle = TRUE, control = list(
#'     burnin = 25, convCrit = "either",
#'     tol0 = 6e-03, tol2 = 6e-03, mcmaxIter = 60))
#' }
bootSE <- function(object, nboot = 100, ci = 0.95, use.mle = TRUE,
                   verbose = FALSE, control = list(), progress = TRUE,
                   ncores = 1, safe.boot = FALSE, ...) {

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

  # Control parameters
  con <- object$control
  nc <- names(con)
  control <- c(control, list(...))
  con[(conArgs <- names(control))] <- control

  if (length(unmatched <- conArgs[!(conArgs %in% nc)]) > 0) {
    warning("Unknown arguments passed to 'control': ", paste(unmatched, collapse = ", "))
  }

  # Use fitted model MLE as initial values?
  if (use.mle) {
    theta.hat <- object$coefficients
    theta.hat <- theta.hat[-which(names(theta.hat) == "haz")]
  } else {
    theta.hat <- NULL
  }

  # Internal function for bootstrap
  if (safe.boot) {
    bootfun <- function() {
      fit.boot <- NULL
      # keep iterating if the model fails with errors
      while (is.null(fit.boot)) {
        # bootstrap sample data
        data.boot <- sampleData(object = object)
        # fit joint model
        fit.boot <- tryCatch(suppressMessages(
          mjoint(formLongFixed = formLongFixed,
                 formLongRandom = formLongRandom,
                 formSurv = formSurv,
                 data = data.boot$longData.boot,
                 survData = data.boot$survData.boot,
                 timeVar = timeVar,
                 inits = theta.hat,
                 verbose = verbose,
                 pfs = FALSE,
                 control = con,
                 ...)),
          error = function(e) NULL)
      }
      return(fit.boot)
    }} else {
      bootfun <- function() {
        # bootstrap sample data
        data.boot <- sampleData(object = object)
        # fit joint model
        fit.boot <- suppressMessages(
          mjoint(formLongFixed = formLongFixed,
                 formLongRandom = formLongRandom,
                 formSurv = formSurv,
                 data = data.boot$longData.boot,
                 survData = data.boot$survData.boot,
                 timeVar = timeVar,
                 inits = theta.hat,
                 verbose = verbose,
                 pfs = FALSE,
                 control = con,
                 ...)
        )
        return(fit.boot)
      }
    }

  if (ncores > 1) {
    ncores.max <- parallel::detectCores()
    if (ncores > ncores.max) {
      ncores <- ncores.max
      warning(paste(
        "ncores is > the maximum number of cores on this machine. Switching to ncores =",
        ncores.max))
    }
    # *** Parallel version ***
    doParallel::registerDoParallel(cores = ncores)
    out <- foreach(b = 1:nboot, .packages = 'joineRML') %dopar% {
      fit.boot <- bootfun()
      return(list("coefs" = fit.boot$coefficient,
                  "conv" = fit.boot$conv))
    }
    registerDoSEQ()

  } else {
    # *** Serial version (incl. progress bar) ***
    out <- list()
    conv.status <- vector(length = nboot)
    if (progress) {
      cat("\n\n") # start progress bar
      pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
    }
    for (b in 1:nboot) {
      fit.boot <- bootfun()
      out[[b]] <- list("coefs" = fit.boot$coefficient,
                       "conv" = fit.boot$conv)
      if (progress) {
        utils::setTxtProgressBar(pb, b) # update progress bar
      }
    }
    if (progress) {
      close(pb) # close progress bar
    }

  }

  conv.status <- sapply(out, function(x) x$conv)

  # Checks
  if (mean(conv.status) <= 0.1) {
    stop("Cannot estimate SEs: fewer than 10% of bootstrap models converged.")
  }
  out <- out[conv.status]

  # Bootstrap estimates
  beta.mat   <- sapply(out, function(u) u$coefs$beta)
  gamma.mat  <- sapply(out, function(u) u$coefs$gamma)
  sigma2.mat <- sapply(out, function(u) u$coefs$sigma2)

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
