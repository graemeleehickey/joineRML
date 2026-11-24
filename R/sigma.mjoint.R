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
#' @details These values are the square roots of the residual variances stored in \code{coef(object)$sigma2}; they are not standard errors of estimated parameters. Names identify the longitudinal outcomes, or use \code{sigma_1}, \code{sigma_2}, and so on when the formula list is unnamed.
#' @importFrom stats sigma
#' @export
sigma.mjoint <- function(object, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  sig <- sqrt(object$coefficients$sigma2)
  sig.names <- names(object$formLongFixed)
  if (is.null(sig.names)) {
    names(sig) <- paste0("sigma_", seq_along(sig))
  } else {
    names(sig) <- sig.names
  }

  sig

}
