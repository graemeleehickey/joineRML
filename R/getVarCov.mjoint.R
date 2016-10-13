#' Extract variance-covariance matrix from an \code{mjoint} object
#'
#' Extract variance-covariance matrix from an \code{mjoint} object.
#'
#' @inheritParams bootSE
#' @param ... further arguments to be passed to or from other methods. They are
#'   ignored in this function.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return A variance-covariance matrix.
#' @export
#'
#' @examples
#' # Fit a joint model with bivariate longitudinal outcomes
#'
#' \dontrun{
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
getVarCov.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- object$coefficients$D
  attr(out, "group.levels") <- object$id

  class(out) <- c("random.effects", "VarCov")
  out

}
