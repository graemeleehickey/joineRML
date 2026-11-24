#' Extract \code{mjoint} residuals
#'
#' @description The residuals at level \emph{i} are obtained by subtracting the fitted
#'   levels at that level from the response vector.
#'
#' @inheritParams confint.mjoint
#' @param level an optional integer giving the level of grouping to be used in
#'   extracting the residuals from object. Level values increase from outermost
#'   to innermost grouping, with level 0 corresponding to the population
#'   residuals and level 1 corresponding to subject-specific residuals. Defaults
#'   to \code{level=0}.
#'
#' @details Subject-specific residuals (\code{level = 1}) require posterior random effects from a converged fit with \code{pfs = TRUE}. Population-level residuals (\code{level = 0}) do not require these post-fit calculations.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link{mjoint}}, \code{\link{fitted.mjoint}}
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return A \code{list} of length \emph{K} with each element a vector of
#'   residuals for the \emph{k}-th longitudinal outcome.
#' @export
residuals.mjoint <- function(object, level = 0, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  dmats <- object$dmats
  if (is.null(dmats)) {
    stop("Need post fit statistics to calculate residuals.
          Re-run the model with 'pfs = TRUE'")
  }

  beta <- object$coefficients$beta

  p <- object$dims$p
  r <- object$dims$r
  K <- object$dims$K
  nik <- object$dmats$l$nik

  Y <- object$dmats$l$yi
  X <- object$dmats$l$Xi
  Z <- object$dmats$l$Zi

  Xbeta <- lapply(X, function(x) {
    x %*% beta
  })
  if (level == 0) {
    ri <- Map(function(y, yhat) {
      y - yhat
    },
    y = Y, yhat = Xbeta)
  } else if (level == 1) {
    .mjoint_require_random_effects(object)
    Eb <- as.data.frame(object$Eb)
    Eb.list <- lapply(rownames(Eb), function(x) Eb[x, ])
    Zb <- Map(function(b, z) z %*% t(b), b = Eb.list, z = Z)
    ri <- Map(function(y, yhat, z) {
      y - yhat - z
    },
    y = Y, yhat = Xbeta, z = Zb)
  } else {
    stop(paste("Unknown level selected:", level))
  }

  resids <- list()
  index <- lapply(nik, function(n) {
    c(0, cumsum(n))
  })
  for (k in 1:K) {
    resids[[k]] <- Map(function(x, i) {
      x[(i[k] + 1):(i[k + 1])]
    },
    x = ri, i = index)
  }
  resids <- lapply(resids, unlist)
  if (!is.null(object$formLongFixed)) {
    names(resids) <- names(object$formLongFixed)
  }

  return(resids)

}
