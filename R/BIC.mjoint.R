#' Extract the Bayesian Information Criterion (BIC) from an \code{mjoint} object
#'
#' @description Extract the Bayesian Information Criterion from an \code{mjoint} object.
#'
#' @param object An object of class \code{mjoint}.
#' @param ... Not used.
#'
#' @author Alessandro Gasparini (\email{ag475@@leicester.ac.uk})
#' @keywords methods
#' @rdname ic
#' @seealso \code{\link[stats]{BIC}} for the generic method description.
#'
#' @return Returns a single number with the BIC value for the \code{object} passed to the function.
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a joint model with bivariate longitudinal outcomes
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
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
#' BIC(fit2)
#' }
BIC.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- summary(object)$BIC
  out
}
