#' Extract the number of subjects
#'
#' @param object A fitted \code{mjoint} object.
#' @param ... Additional arguments (not used).
#' @return The number of subjects used to fit the joint model. Longitudinal observation counts are reported separately by \code{glance()}.
#' @details This agrees with the \code{nobs} attribute of \code{logLik(object)} and the sample size used in BIC.
#' @author Graeme L. Hickey
#' @export
nobs.mjoint <- function(object, ...) {
  object$dims$n
}

# Validate the outcome selection consistently across extraction methods.
.mjoint_outcome <- function(object, k) {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
      !is.finite(k) || k != floor(k) || k < 1 || k > object$dims$K) {
    stop("'k' must be an integer between 1 and the number of longitudinal outcomes.")
  }
  as.integer(k)
}

.mjoint_check_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    labels <- names(dots)
    if (is.null(labels)) labels <- rep("", length(dots))
    labels[is.na(labels) | labels == ""] <- "<unnamed>"
    stop("Unused argument(s): ", paste(labels, collapse = ", "), ".")
  }
}

.mjoint_outcome_names <- function(object) {
  labels <- names(object$formLongFixed)
  if (is.null(labels)) labels <- rep("", object$dims$K)
  for (k in seq_len(object$dims$K)) {
    if (is.na(labels[k]) || labels[k] == "") {
      labels[k] <- paste(deparse(object$formLongFixed[[k]][[2]]), collapse = "")
    }
  }
  labels
}

#' Extract model terms, frames, and design matrices
#'
#' @param x,object,formula A fitted \code{mjoint} object.
#' @param process Select \code{"Longitudinal"} (the default) or \code{"Event"}.
#' @param k The longitudinal outcome number, between 1 and the number of outcomes. Only used for \code{process = "Longitudinal"}.
#' @param ... Additional arguments. These methods extract the fitted model only; supplying other arguments, including new data, is an error.
#' @details The longitudinal methods describe the fixed-effects model for outcome \code{k}. The event methods describe the baseline covariates in the survival formula. Its design matrix has no intercept or latent association columns. Rows follow the stored longitudinal subject-time order or survival subject order.
#'
#' Model frames contain the evaluated formula variables, including transformed responses and covariates. Fitted terms and contrasts are retained, so changing the original data or R's contrast options does not change extraction. For objects fitted before version 0.4.9, frames are reconstructed from stored data and terms when available.
#'
#' Column names of a longitudinal design matrix match \code{coef(object, process = "Longitudinal", k = k)}. Event design columns match the baseline-covariate entries of \code{coef(object, process = "Event")}; the latter also includes association coefficients.
#' @return \code{terms()} returns a \code{terms} object; \code{model.frame()} returns a data frame; \code{model.matrix()} returns a numeric matrix.
#' @author Graeme L. Hickey
#' @name model-extraction
NULL

#' @rdname model-extraction
#' @export
model.frame.mjoint <- function(formula, process = c("Longitudinal", "Event"),
                               k = 1, ...) {
  .mjoint_check_dots(...)
  object <- formula
  process <- match.arg(process)
  if (process == "Longitudinal") {
    k <- .mjoint_outcome(object, k)
    frame <- object$model.frames$longitudinal[[k]]
    if (is.null(frame)) {
      frame <- model.frame(object$lfit[[k]]$terms, object$data[[k]],
                           na.action = na.fail)
      attr(frame, "contrasts") <- object$lfit[[k]]$contrasts
    }
  } else {
    frame <- object$model.frames$event
    if (is.null(frame)) {
      frame <- model.frame(object$sfit$terms, object$survData,
                           na.action = na.fail, xlev = object$sfit$xlevels)
    }
  }
  frame
}

#' @rdname model-extraction
#' @export
terms.mjoint <- function(x, process = c("Longitudinal", "Event"), k = 1, ...) {
  attr(model.frame(x, process = match.arg(process), k = k, ...), "terms")
}

#' @rdname model-extraction
#' @export
model.matrix.mjoint <- function(object, process = c("Longitudinal", "Event"),
                                k = 1, ...) {
  .mjoint_check_dots(...)
  process <- match.arg(process)
  if (process == "Event") return(model.matrix(object$sfit))
  frame <- model.frame(object, process = process, k = k)
  contrasts <- attr(frame, "contrasts")
  contrasts <- contrasts[intersect(names(contrasts), names(frame))]
  model.matrix(attr(frame, "terms"), frame, contrasts.arg = contrasts)
}

#' Extract joint-model coefficients
#'
#' @param object A fitted \code{mjoint} object.
#' @param complete Logical; retained for compatibility with \code{coef.default()}.
#' @param process \code{NULL} (the default) returns the original coefficient list. Select \code{"Longitudinal"} or \code{"Event"} to return a named vector of regression coefficients.
#' @param k With \code{process = "Longitudinal"}, select one outcome. The default \code{NULL} returns coefficients for all longitudinal outcomes, retaining their numbered suffixes. Selecting an outcome removes those suffixes so names match its design matrix.
#' @param ... Additional arguments (not supported).
#' @details The default return value is unchanged. Event coefficients include both baseline covariates and latent association coefficients. The same \code{process} and \code{k} selections can be supplied to \code{vcov()} to obtain the corresponding covariance matrix in the same order.
#' @return A list of all estimated parameters, or a named numeric vector for a selected component.
#' @author Graeme L. Hickey
#' @export
coef.mjoint <- function(object, complete = TRUE, ..., process = NULL, k = NULL) {
  .mjoint_check_dots(...)
  if (is.null(process)) {
    if (!is.null(k)) stop("Specify 'process = \"Longitudinal\"' when selecting 'k'.")
    out <- object$coefficients
  } else {
    selection <- .mjoint_coef_selection(object, process, k)
    out <- c(object$coefficients$beta, object$coefficients$gamma)[selection$index]
    names(out) <- selection$names
  }
  if (complete) out else out[!is.na(out)]
}

.mjoint_coef_selection <- function(object, process, k = NULL) {
  process <- match.arg(process, c("Longitudinal", "Event"))
  nb <- sum(object$dims$p)
  if (process == "Event") {
    if (!is.null(k)) stop("'k' is only used for longitudinal coefficients.")
    index <- nb + seq_along(object$coefficients$gamma)
    labels <- names(object$coefficients$gamma)
  } else if (is.null(k)) {
    index <- seq_len(nb)
    labels <- names(object$coefficients$beta)
  } else {
    k <- .mjoint_outcome(object, k)
    index <- cumsum(c(0, object$dims$p))[k] + seq_len(object$dims$p[k])
    labels <- sub(paste0("_", k, "$"), "", names(object$coefficients$beta)[index])
  }
  list(index = index, names = labels)
}
