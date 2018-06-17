#' @title Extract information criteria from an \code{mjoint} object
#'
#' @description Extract the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC) from an \code{mjoint} object.
#'
#' @param object An object of class \code{mjoint}.
#' @param ... Not used.
#' @param k The penalty applied to the AIC. Defaults to 2.
#'
#' @author Alessandro Gasparini (\email{ag475@@leicester.ac.uk})
#' @keywords methods
#' @rdname ic
#' @seealso \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}} for the generic methods description.
#'
#' @return Returns a single number with the AIC or the BIC value for the \code{object} passed to the function.
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a joint model with bivariate longitudinal outcomes
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) &
#'                    !is.na(heart.valve$log.lvmi), ]
#' fit2 <- mjoint(
#'     formLongFixed = list("grad" = log.grad ~ time + sex + hs,
#'                          "lvmi" = log.lvmi ~ time + sex),
#'     formLongRandom = list("grad" = ~ 1 | num,
#'                           "lvmi" = ~ time | num),
#'     formSurv = Surv(fuyrs, status) ~ age,
#'     data = list(hvd, hvd),
#'     inits = list("gamma" = c(0.11, 1.51, 0.80)),
#'     timeVar = "time")
#'
#' AIC(fit2)
#' BIC(fit2)
#' }
AIC.mjoint <- function(object, ..., k = 2) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- k * nrow(vcov(object)) - 2 * as.numeric(logLik(object))
  out
}

#' @rdname ic
#' @note The penalty value for the BIC is the natural logarithm of the number of individuals.
#' @export
BIC.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- AIC(object, k = log(object$dims$n))
  out
}
