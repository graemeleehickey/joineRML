#' Extract variance-covariance matrix of random effects from an \code{mjoint}
#' object
#'
#' @description Extract variance-covariance matrix of random effects from an
#'   \code{mjoint} object.
#'
#' @param obj an object inheriting from class \code{mjoint} for a joint model of
#'   time-to-event and multivariate longitudinal data.
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[nlme]{getVarCov}} for the generic method description.
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return A variance-covariance matrix.
#' @importFrom nlme getVarCov
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
#' getVarCov(fit2)
#' }
getVarCov.mjoint <- function(obj, ...) {

  if (!inherits(obj, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- obj$coefficients$D
  attr(out, "group.levels") <- obj$id

  class(out) <- c("random.effects", "VarCov")
  out

}
