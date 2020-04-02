#' Extract residual standard deviation(s) from an \code{mjoint} object
#'
#' @description Extract residual standard deviation(s) from an \code{mjoint}
#'   object.
#'
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[lme4]{sigma}} in the \strong{lme4} package.
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return a number (standard deviation) if \eqn{K = 1} (univariate model), or a
#'   vector if \eqn{K>1} (multivariate model).
#' @importFrom stats sigma
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
