#' Extract log-likelihood from an \code{mjoint} object
#'
#' @description Extract log-likelihood from an \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[stats]{logLik}} for the generic method description.
#'
#' @references
#'
#' Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
#' measurements and event time data. \emph{Biostatistics.} 2000; \strong{1(4)}:
#' 465-480.
#'
#' @return Returns an object of class \code{logLik}. This is a number with two
#'   attributes: \code{df} (degrees of freedom), giving the number of parameters
#'   in the model, and \code{nobs}, the number of observations used in
#'   estimation.
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
#'     timeVar = "time",
#'     verbose = TRUE)
#'
#' logLik(fit2)
#' }
logLik.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- object$log.lik
  attr(out, "df") <- with(object$dims, sum(p) + q + 2*K + (sum(r) * (sum(r) + 1)) / 2)
  attr(out, "nobs") <- sum(object$dims$nk)

  class(out) <- "logLik"
  out

}
