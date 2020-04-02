#' Extract model formulae from an \code{mjoint} object
#'
#' @description Extract model formulae from an \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#' @param x an object inheriting from class \code{mjoint} for a joint model of
#'   time-to-event and multivariate longitudinal data.
#' @param process character string: if \code{process='Longitudinal'} a fixed
#'   effects formula from the (multivariate) longitudinal sub-model is returned
#'   for the \code{k}-th outcome. Else, if \code{process='Event'}, the
#'   time-to-event model formula is returned.
#' @param k integer: a number between 1 and \emph{K} (the total number of
#'   longitudinal outcomes) that specifies the longitudinal outcome of interest.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[stats]{formula}} for the generic method description, and
#'   \code{\link{ranef.mjoint}}.
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
#' measured with error. \emph{Biometrics.} 1997; \strong{53(1)}: 330-339.
#'
#' @return An object of class "formula" which contains a symbolic model formula
#'   for the separate sub-model fixed effect terms only.
#' @export
formula.mjoint <- function(x, process = c("Longitudinal", "Event"), k = 1, ...) {

  if (!inherits(x, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  K <- x$dims$K
  process <- match.arg(process)

  formOut <- NULL

  if (process == "Longitudinal") {
    if (is.null(k) || is.na(k)) {
      stop("Must specify a longitudinal outcome.")
    }
    if (k > K) {
      stop("Incompatible with dimensions of the joint model.")
    }
    formOut <- x$formLongFixed[[k]]
  } else {
    formOut <- x$formSurv
  }

  class(formOut) <- "formula"
  formOut

}
