#' Extract \code{mjoint} fitted values
#'
#' @description The fitted values at level \emph{i} are obtained by adding
#'   together the population fitted values (based only on the fixed effects
#'   estimates) and the estimated contributions of the random effects.
#'
#' @inheritParams confint.mjoint
#' @param level an optional integer giving the level of grouping to be used in
#'   extracting the fitted values from object. Level values increase from outermost
#'   to innermost grouping, with level 0 corresponding to the population
#'   fitted values and level 1 corresponding to subject-specific fitted values Defaults
#'   to level 0.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link{mjoint}}, \code{\link{residuals.mjoint}}
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return A \code{list} of length \emph{K} with each element a vector of
#'   fitted values for the \emph{k}-th longitudinal outcome.
#' @export
fitted.mjoint <- function(object, level = 0, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  dmats <- object$dmats
  if (is.null(dmats)) {
    stop("Need post fit statistics to calculate fitted values
         Re-run the model with 'pfs = TRUE'")
  }

  beta <- object$coefficients$beta
  Eb <- as.data.frame(object$Eb)

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
  Eb.list <- lapply(rownames(Eb), function(x) {
    Eb[x, ]
  })
  Zb <- mapply(function(b, z) {
    z %*% t(b)
  },
  b = Eb.list, z = Z)

  if (level == 0) {
    ri <- Xbeta
  } else if (level == 1) {
    ri <- mapply(function(y, yhat, z) {
      yhat + z
    },
    yhat = Xbeta, z = Zb)
  } else {
    stop(paste("Unknown level selected:", level))
  }

  fvals <- list()
  index <- lapply(nik, function(n) {
    c(0, cumsum(n))
  })
  for (k in 1:K) {
    fvals[[k]] <- mapply(function(x, i) {
      x[(i[k] + 1):(i[k + 1])]
    },
    x = ri, i = index)
  }
  fvals <- lapply(fvals, unlist)
  if (!is.null(object$formLongFixed)) {
    names(fvals) <- names(object$formLongFixed)
  }

  return(fvals)

  }
