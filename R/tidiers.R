#' Tidying methods for joint models for time-to-event data and multivariate
#' longitudinal data
#'
#' These methods tidy the coefficients of joint models for time-to-event data
#' and multivariate longitudinal data of the \code{mjoint} class from the
#' \code{joineRML} package.
#'
#' @param x An object of class \code{mjoint}.
#'
#' @return All tidying methods return a \code{data.frame} without rownames. The
#'   structure depends on the method chosen.
#'
#' @name mjoint_tidiers
#'
#' @examples
#' \dontrun{
#' # Fit a joint model with bivariate longitudinal outcomes
#' library(joineRML)
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) &
#'   !is.na(heart.valve$log.lvmi) &
#'   heart.valve$num <= 50, ]
#' fit <- mjoint(
#'   formLongFixed = list(
#'     "grad" = log.grad ~ time + sex + hs,
#'     "lvmi" = log.lvmi ~ time + sex
#'   ),
#'   formLongRandom = list(
#'     "grad" = ~ 1 | num,
#'     "lvmi" = ~ time | num
#'   ),
#'   formSurv = Surv(fuyrs, status) ~ age,
#'   data = hvd,
#'   inits = list("gamma" = c(0.11, 1.51, 0.80)),
#'   timeVar = "time"
#' )
#'
#' # Extract the survival fixed effects
#' tidy(fit)
#'
#' # Extract the longitudinal fixed effects
#' tidy(fit, component = "longitudinal")
#'
#' # Extract the survival fixed effects with confidence intervals
#' tidy(fit, conf.int = TRUE)
#'
#' # Extract the survival fixed effects with confidence intervals based on 
#' # bootstrapped standard errors
#' bSE <- bootSE(fit, nboot = 5, safe.boot = TRUE)
#' tidy(fit, bootSE = bSE, conf.int = TRUE)
#'
#' # Augment original data with fitted longitudinal values and residuals
#' hvd2 <- augment(fit)
#'
#' # Extract model statistics
#' glance(fit)
#' }
#'
#' @rdname mjoint_tidiers
#'
#' @param component Select \code{"survival"} (the default), \code{"longitudinal"}, or \code{"all"}.
#'
#' @param bootSE An object of class \code{bootSE} for the corresponding model.
#'   If \code{bootSE = NULL} (the default), the function will use approximate
#'   standard error estimates calculated from the empirical information matrix.
#'
#' @param conf.int Include (1 - \code{conf.level})\% confidence intervals?
#'   Defaults to \code{FALSE}.
#'
#' @param conf.level The confidence level required.
#' @param boot_se An alias for \code{bootSE}, matching broom's argument name. Supply only one of these arguments.
#' @param exponentiate Logical; exponentiate survival estimates and confidence limits. Standard errors and test statistics remain on the log-hazard scale. Only supported for \code{component = "survival"}. Exponentiated association coefficients refer to a one-unit change in the associated latent process, not a one-unit change in the observed marker.
#'
#' @author Alessandro Gasparini (\email{alessandro.gasparini@@ki.se})
#' 
#' @return \code{tidy} returns one row for each estimated fixed effect depending
#'   on the \code{component} parameter. It contains the following  columns:
#'   \item{term}{The term being estimated} \item{estimate}{Estimated value}
#'   \item{std.error}{Standard error} \item{statistic}{Z-statistic}
#'   \item{p.value}{P-value computed from Z-statistic} \item{conf.low}{The lower
#'   bound of a confidence interval on \code{estimate}, if required}
#'   \item{conf.high}{The upper bound of a confidence interval on
#'   \code{estimate}, if required}. Additional columns identify the \code{component}, \code{outcome}, and \code{effect} (\code{"fixed"} or \code{"association"}). Outcome labels use the names of \code{formLongFixed}, falling back to the response expressions. Survival covariates are labelled \code{"event"}; association coefficients are labelled by their corresponding longitudinal outcome. Existing term names are retained.
#'
#' @export
#' @exportS3Method broom::tidy
tidy.mjoint <- function(x, 
                        component = "survival", 
                        bootSE = NULL, 
                        conf.int = FALSE, 
                        conf.level = 0.95, 
                        ..., boot_se = NULL, exponentiate = FALSE) {
  
  .mjoint_check_dots(...)
  component <- match.arg(component, c("survival", "longitudinal", "all"))
  for (flag in c("conf.int", "exponentiate")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("'", flag, "' must be TRUE or FALSE.")
    }
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a number strictly between 0 and 1.")
  }
  if (exponentiate && component != "survival") {
    stop("Exponentiation is only supported for 'component = \"survival\"'.")
  }
  if (!is.null(bootSE) && !is.null(boot_se)) {
    stop("Supply only one of 'bootSE' and 'boot_se'.")
  }
  if (!is.null(boot_se)) bootSE <- boot_se
  
  if (!is.null(bootSE)) {
    if (!inherits(x = bootSE, what = "bootSE")) stop("'bootSE' object not of class 'bootSE'")
  }

  # make summary object
  if (is.null(bootSE)) {
    smr <- summary(x)
  } else {
    smr <- summary(x, bootSE = bootSE)
  }

  # extract appropriate component
  labels <- .mjoint_outcome_names(x)
  metadata <- data.frame(
    component = c(rep("longitudinal", sum(x$dims$p)),
                  rep("survival", x$dims$q + x$dims$K)),
    outcome = c(rep(labels, x$dims$p), rep("event", x$dims$q), labels),
    effect = c(rep("fixed", sum(x$dims$p) + x$dims$q),
               rep("association", x$dims$K))
  )
  if (component == "survival") {
    out <- data.frame(smr$coefs.surv)
    term <- rownames(smr$coefs.surv)
  } else if (component == "longitudinal") {
    out <- data.frame(smr$coefs.long)
    term <- rownames(smr$coefs.long)
  } else {
    out <- data.frame(rbind(smr$coefs.long, smr$coefs.surv))
    term <- c(rownames(smr$coefs.long), rownames(smr$coefs.surv))
  }
  if (component != "all") metadata <- metadata[metadata$component == component, ]

  # fix names
  names(out) <- c("estimate", "std.error", "statistic", "p.value")

  # turn rownames into a 'term' column
  out$term <- term
  rownames(out) <- NULL
  out <- out[, c(5, 1:4)]

  # make confidence intervals (if required)
  if (conf.int) {
    cv <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - cv * out$std.error
    out$conf.high <- out$estimate + cv * out$std.error
  }
  if (exponentiate) {
    columns <- intersect(c("estimate", "conf.low", "conf.high"), names(out))
    out[columns] <- lapply(out[columns], exp)
  }
  out <- cbind(out, metadata)

  # turn out into a tibble object
  out <- tibble::as_tibble(out)

  # return tidy object
  return(out)
  
}


#' @rdname mjoint_tidiers
#'
#' @param data Original data as a data frame or list of identical data frames. Extracted from \code{x} by default. Supplied records are matched by subject ID and measurement time and returned in the supplied order; these pairs must be unique when reordering is required. All original columns must be present and unchanged. Different observation schedules across outcomes cannot be combined into a single augmented data frame.
#'   
#' @author Alessandro Gasparini (\email{alessandro.gasparini@@ki.se})
#' 
#' @return \code{augment} returns one row for each original observation, with
#'   columns (each prepended by a .) added. Included are the columns:
#'   \item{.fitted_j_0}{population-level fitted values for the j-th longitudinal
#'   process} \item{.fitted_j_1}{individuals-level fitted values for the j-th
#'   longitudinal process} \item{.resid_j_0}{population-level residuals for the
#'   j-th longitudinal process} \item{.resid_j_1}{individual-level residuals for
#'   the j-th longitudinal process} See \code{\link[joineRML]{fitted.mjoint}}
#'   and \code{\link[joineRML]{residuals.mjoint}} for more information on the
#'   difference between population-level and individual-level fitted values and
#'   residuals.
#'
#' @note Augmented column names use the longitudinal outcome labels. Repeated labels are made unique by adding numbered suffixes.
#'
#' @export
#' @exportS3Method broom::augment
augment.mjoint <- function(x, data = x$data, ...) {
  .mjoint_check_dots(...)
 
  # checks on 'data'
  if (is.null(data) || is.null(x$data)) {
    stop("The original data are unavailable in 'x'; observations cannot be matched.")
  }

  if (!is.data.frame(data)) {
    if (!is.list(data) || !length(data) ||
        !all(vapply(data, is.data.frame, logical(1)))) {
      stop("'data' must be a data frame or a list of data frames.")
    }
    if (!all(vapply(data, identical, logical(1), data[[1]]))) {
      stop("List of 'data' extracted from 'x' does not include identical data frames.")
    }
    data <- data[[1]]
  }
  if (!all(vapply(x$data, identical, logical(1), x$data[[1]]))) {
    stop("Cannot augment outcomes with different data frames or observation schedules.")
  }
  original <- x$data[[1]]
  if (nrow(data) != nrow(original) || !all(names(original) %in% names(data))) {
    stop("'data' must contain all original observations and columns.")
  }
  same <- function(a, b) identical(lapply(a, as.character), lapply(b, as.character))
  row.index <- seq_len(nrow(original))
  if (!same(data[names(original)], original)) {
    keys <- c(x$id, x$timeVar[1])
    left <- data[keys]
    right <- original[keys]
    left[[x$id]] <- as.character(left[[x$id]])
    right[[x$id]] <- as.character(right[[x$id]])
    if (anyDuplicated(left) || anyDuplicated(right)) {
      stop("Cannot reorder 'data' with duplicate subject-time pairs; use augment(x).")
    }
    left$.input_row <- seq_len(nrow(left))
    right$.fitted_row <- seq_len(nrow(right))
    matched <- merge(left, right, by = keys, all.x = TRUE, sort = FALSE)
    row.index <- matched$.fitted_row[order(matched$.input_row)]
    if (anyNA(row.index) || !same(data[names(original)], original[row.index, , drop = FALSE])) {
      stop("'data' must contain unchanged original observations; new data are not supported.")
    }
  }

  # longitudinal fitted values
  labels <- make.unique(.mjoint_outcome_names(x), sep = "_")
  fit0 <- fitted(x, level = 0)
  names(fit0) <- paste0(".fitted_", labels, "_0")
  fit0 <- do.call(cbind.data.frame, fit0)
  
  fit1 <- fitted(x, level = 1)
  names(fit1) <- paste0(".fitted_", labels, "_1")
  fit1 <- do.call(cbind.data.frame, fit1)

  # longitudinal residuals
  res0 <- residuals(x, level = 0)
  names(res0) <- paste0(".resid_", labels, "_0")
  res0 <- do.call(cbind.data.frame, res0)
  res1 <- residuals(x, level = 1)
  names(res1) <- paste0(".resid_", labels, "_1")
  res1 <- do.call(cbind.data.frame, res1)

  # return augmented 'data'
  additions <- cbind(fit0, fit1, res0, res1)
  out <- cbind(data, additions[row.index, , drop = FALSE])
  out <- tibble::as_tibble(out)
  return(out)
  
}


#' @rdname mjoint_tidiers
#'
#' @param ... Additional arguments. Unsupported arguments cause an error, including \code{newdata} for \code{augment()}.
#'
#' @author Alessandro Gasparini (\email{alessandro.gasparini@@ki.se})
#' 
#' @return \code{glance} returns one row with the columns \item{sigma_j}{the
#'   estimated residual standard deviation for the j-th longitudinal
#'   process} \item{AIC}{the Akaike Information Criterion} \item{BIC}{the
#'   Bayesian Information Criterion} \item{logLik}{the data's log-likelihood
#'   under the model}. Additional columns are \code{nobs} (subjects), \code{nevent} (observed events), \code{npar} (the parameter count used by AIC and BIC, excluding baseline hazard increments), \code{converged}, and \code{nobs_longitudinal_j} (measurement counts for outcome \emph{j}). This summary is available without post-fit standard errors.
#'
#' @export
#' @exportS3Method broom::glance
glance.mjoint <- function(x, ...) {
  .mjoint_check_dots(...)
  ll <- logLik(x)
  out <- as.data.frame(as.list(setNames(sigma(x),
                                       paste0("sigma_", seq_len(x$dims$K)))))
  out$AIC <- AIC(x)
  out$BIC <- BIC(x)
  out$logLik <- as.numeric(ll)
  out$nobs <- nobs(x)
  out$nevent <- x$sfit$nevent
  out$npar <- attr(ll, "df")
  out$converged <- isTRUE(x$conv)
  out[paste0("nobs_longitudinal_", seq_len(x$dims$K))] <- as.list(x$dims$nk)
  rownames(out) <- NULL
  out <- tibble::as_tibble(out)
  return(out)
  
}
