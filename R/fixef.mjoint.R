#' Extract fixed effects estimates from an \code{mjoint} object
#'
#' Extract fixed effects estimates from an \code{mjoint} object.
#'
#' @param object an object inheriting from class \code{mjoint} for a joint model
#'   of time-to-event and multivariate longitudinal data.
#' @param process character string: if \code{process='Longitudinal'} the fixed
#'   effects coefficients from (multivariate) longitudinal sub-model are
#'   returned. Else, if \code{process='Event'}, the coefficients from the
#'   time-to-event sub-model are returned.
#' @param ... additional arguments; currently none are used.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods
#' @seealso \code{\link[nlme]{fixef}} for the generic method description, and
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
#' @return A named vector of length equal to the number of sub-model coefficients estimated.
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
#' fixef(fit2, process = "Longitudinal")
#' fixef(fit2, process = "Event")
#' }
fixef.mjoint <- function(object, process = c("Longitudinal", "Event"), ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  process <- match.arg(process)
  if (process == "Longitudinal") {
    object$coefficients$beta
  } else {
    object$coefficients$gamma
  }

}