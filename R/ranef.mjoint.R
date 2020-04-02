#' Extract random effects estimates from an \code{mjoint} object
#'
#' @description Extract random effects estimates from an \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#' @param postVar logical: if \code{TRUE} the variance of the posterior
#'   distribution is also returned.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[nlme]{ranef}} for the generic method description, and
#'   \code{\link{fixef.mjoint}}. To plot \code{ranef.mjoint} objects, see
#'   \code{\link{plot.ranef.mjoint}}.
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
#' measured with error. \emph{Biometrics.} 1997; \strong{53(1)}: 330-339.
#'
#' @return A \code{data.frame} (also of class \code{ranef.mjoint}) with rows
#'   denoting the individuals and columns the random effects (e.g., intercepts,
#'   slopes, etc.). If \code{postVar=TRUE}, the numeric matrix has an extra
#'   attribute, \code{postVar}.
#' @importFrom lme4 ranef
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
#' ranef(fit2)
#' }
ranef.mjoint <- function(object, postVar = FALSE, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- object$Eb
  out <- as.data.frame(out)

  if (postVar) {
    attr(out, "postVar") <- object$Vb
  }

  class(out) <- c("ranef.mjoint", "data.frame")
  return(out)

}
