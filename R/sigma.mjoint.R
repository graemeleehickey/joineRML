#' Extract residual standard deviation(s) from an \code{mjoint} object
#'
#' Extract residual standard deviation(s) from an \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods
#' @seealso \code{\link[lme4]{sigma}} in the \strong{lme4} package.
#'
#' @rawNamespace if (getRversion() >= '3.3.0') importFrom(stats,sigma) else
#'   importFrom(lme4,sigma)
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return a number (standard deviation) if \eqn{K=1} (univariate model), or a
#'   vector if \eqn{K>1} (multivariate model).
#' @export
sigma.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  sig <- sqrt(object$coef$sigma2)
  sig.names <- names(object$formLongFixed)
  if (is.null(sig.names)) {
    names(sig) <- paste0("sigma_", 1:length(sig))
  } else {
    names(sig) <- sig.names
  }

  sig

}
