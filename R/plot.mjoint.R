#' Plot diagnostics from an \code{mjoint} object
#'
#' @description Plot diagnostics from an \code{mjoint} object.
#'
#' @param x an object inheriting from class \code{mjoint} for a joint model of
#'   time-to-event and multivariate longitudinal data.
#' @param type currently the only option is \code{type='convergence'} for
#'   graphical examination of convergence over MCEM iteration.
#' @param ... other parameters passed to \code{\link{plotConvergence}}.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods dplot
#' @seealso \code{\link[graphics]{plot.default}}, \code{\link[graphics]{par}},
#'   \code{\link[graphics]{abline}}.
#'
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
#' plot(fit1, param = "beta")  # LMM fixed effect parameters
#' plot(fit1, param = "gamma") # event model parameters
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
#'     control = list(burnin = 50),
#'     verbose = TRUE)
#'
#' plot(fit2, type = "convergence", params = "gamma")
#' }
plot.mjoint <- function(x, type = "convergence", ...) {

  if (!inherits(x, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  if (type == "convergence") {
    plotConvergence(x, ...)
  }

}
