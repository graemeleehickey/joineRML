#' Extract an approximate variance-covariance matrix of estimated parameters
#' from an \code{mjoint} object
#'
#' @description Returns the variance-covariance matrix of the main parameters of
#'   a fitted \code{mjoint} model object.
#'
#' @inheritParams confint.mjoint
#' @param correlation logical: if \code{TRUE} returns the correlation matrix,
#'   otherwise returns the variance-covariance matrix (default).
#'
#' @details This is a generic function that extracts the variance-covariance
#'   matrix of parameters from an \code{mjoint} model fit. It is based on a
#'   profile likelihood, so no estimates are given for the baseline hazard
#'   function, which is generally considered a nuisance parameter. It is based
#'   on the empirical information matrix (see Lin et al. 2002, and McLachlan
#'   and Krishnan 2008 for details), so is only approximate.
#'
#' @note This function is not to be confused with \code{\link{getVarCov}}, which
#'   returns the extracted variance-covariance matrix for the random effects
#'   distribution.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link[stats]{vcov}} for the generic method description, and
#'   \code{\link[stats]{cov2cor}} for details of efficient scaling of a
#'   covariance matrix into the corresponding correlation matrix.
#'
#' @references
#'
#' Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the joint
#' analysis of time-to-event and multiple longitudinal variables. \emph{Stat
#' Med.} 2002; \strong{21}: 2369-2382.
#'
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions}. Second
#' Edition. Wiley-Interscience; 2008.
#'
#' @import stats
#' @importFrom MASS ginv
#'
#' @return A variance-covariance matrix.
#' @export
#'
#' @examples
#' # Fit a classical univariate joint model with a single longitudinal outcome
#' # and a single time-to-event outcome
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#'
#' set.seed(1)
#' fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
#'     formLongRandom = ~ time | num,
#'     formSurv = Surv(fuyrs, status) ~ age,
#'     data = hvd,
#'     timeVar = "time",
#'     control = list(nMCscale = 2, burnin = 5)) # controls for illustration only
#'
#' vcov(fit1)
#'
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
#' vcov(fit2)
#' }
vcov.mjoint <- function(object, correlation = FALSE, ...) {

  # Have used code by Dimitris Rizopoulos here after having some
  # issues inverting the matrix.
  # URL: https://github.com/drizopoulos/JM/blob/master/R/vcov.jointModel.R

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  out <- try(qr.solve(object$Hessian, tol = 1e-12), silent = TRUE)

  if (!inherits(out, "try-error")) {
    vmat <- structure(out, dimnames = dimnames(object$Hessian))
  } else {
    vmat <- structure(MASS::ginv(object$Hessian), dimnames = dimnames(object$Hessian))
  }
  vmat <- (vmat + t(vmat)) / 2

  if (!correlation) {
    vmat
  } else {
    stats::cov2cor(vmat)
  }

}
