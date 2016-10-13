#' Plot diagnostics from an \code{mjoint} object
#'
#' Plot diagnostics from an \code{mjoint} object.
#'
#' @inheritParams bootSE
#' @param type currently the only option is \code{type='convergence'} for
#'   graphical examination of convergence over MCEM iteration.
#' @param ... other parameters passed to \code{\link{plotConvergence}}.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods dplot
#'
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
#' plot(fit2, type = "convergence", params = "gamma")
#' }
plot.mjoint <- function(object, type = "convergence", ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  if (type == "convergence") {
    plotConvergence(object, ...)
  }

}
