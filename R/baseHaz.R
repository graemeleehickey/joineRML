#' The baseline hazard estimate of an \code{mjoint} object
#'
#' @description This function returns the (baseline) hazard increment from a
#'   fitted \code{mjoint} object. In addition, it can report either the
#'   \emph{uncentered} or the more ubiquitous \emph{centered} version.
#'
#' @inheritParams confint.mjoint
#' @param centered logical: should the baseline hazard be for the mean-centered
#'   covariates model or not? Default is \code{centered=TRUE}. See
#'   \strong{Details}.
#'
#' @details When covariates are included in the time-to-event sub-model,
#'   \code{\link{mjoint}} automatically centers them about their respective
#'   means. This also applies to non-continuous covariates, which are first
#'   coded using a dummy-transformation for the design matrix and subsequently
#'   centered. The reason for the mean-centering is to improve numerial
#'   statbility, as the survival function involves exponential terms. Extracting
#'   the baseline hazard increments from \code{\link{mjoint.object}} returns the
#'   Breslow hazard estimate (Lin, 2007) that corresponds to this mean-centered
#'   model. This is the same as is done in the R \code{survival} package when
#'   using \code{\link[survival]{coxph.detail}} (Therneau and Grambsh, 2000). If
#'   the user wants to access the baseline hazard estimate for the model in
#'   which no mean-centering is applied, then they can use this function, which
#'   scales the mean-centered baseline hazard by
#'
#'   \deqn{\exp\{-\bar{w}^\top \gamma_v\},}
#'
#'   where \eqn{\bar{w}} is a vector of the means from the time-to-event
#'   sub-model design matrix.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods survival
#' @seealso \code{\link{mjoint}} and \code{\link[stats]{coef}}.
#'
#' @references
#'
#' Therneau TM, Grambsch PM. \emph{Modeling Survival Data: Extending the Cox
#' Model.} New Jersey: Springer-Verlag; 2000.
#'
#' Lin DY. On the Breslow estimator. \emph{Lifetime Data Anal.} 2007;
#' \strong{13(4)}: 471-480.
#'
#' @return A \code{data.frame} with 2 columns: the unique failure times and the
#'   estimate baseline hazard.
#' @export
#'
#' @examples
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
#' baseHaz(fit2, centered = FALSE)
#' }
baseHaz <- function(object, centered = TRUE) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  times <- object$dmats$t$tj
  haz <- object$coefficients$haz
  q <- object$dims$q

  if ((q == 0) && centered) {
    warning("No covariates in model to centre.\n")
  }

  if ((q > 0) && !centered) {
    xcenter <- object$dmats$t$xcenter
    gamma.v <- object$coefficients$gamma[1:q]
    haz <- haz * exp(-sum(xcenter * gamma.v))
  }

  out <- data.frame("time" = times, "haz" = haz)
  return(out)

}
