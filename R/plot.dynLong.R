#' Plot a \code{dynLong} object
#'
#' @description Plots the conditional longitudinal expectations for a
#'   \emph{new} subject calculated using the \code{\link{dynLong}} function.
#'
#' @param x an object of class \code{dynLong} calculated by the
#'   \code{\link{dynLong}} function.
#' @inheritParams plot.dynSurv
#' @param ylab a character vector of the titles for the \emph{K} longitudinal
#'   outcomes y-axes: see \code{\link[graphics]{title}}.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords hplot
#' @seealso \code{\link{dynLong}}
#'
#' @references
#'
#' Rizopoulos D. Dynamic predictions and prospective accuracy in joint models
#' for longitudinal and time-to-event data. \emph{Biometrics}. 2011;
#' \strong{67}: 819–829.
#'
#' @return A dynamic prediction plot.
#' @import graphics
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
#' hvd2 <- droplevels(hvd[hvd$num == 1, ])
#' out <- dynLong(fit2, hvd2)
#' plot(out, main = "Patient 1")
#' }
plot.dynLong <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                         grid = TRUE, ...) {

  if (!inherits(x, "dynLong")) {
    stop("Use only with 'dynLong' model objects.\n")
  }

  # Extract information we need
  fit <- x$fit
  pred <- x$pred
  newdata <- x$newdata
  data.t <- x$data.t

  # Set-up the plotting grid
  old.par <- par(no.readonly = TRUE)
  K <- fit$dims$K
  par(mfrow = c(K, 1), mar = c(0, 4.5, 0, 2), oma = c(4, 0, 3, 0))

  if (!is.null(ylab)) {
    if (length(ylab) != K) {
      stop("Number of longitudinal axes labels does not match number of outcomes.\n")
    }
  }

  ylimfun <- function(k) {
    y.min <- min(data.t$yk[[k]], pred[[k]][, 2])
    y.max <- max(data.t$yk[[k]], pred[[k]][, 2])
    return(c(y.min, y.max))
  }

  # Plot longitudinal markers (extrapolated)
  for (k in 1:K) {
    plot(pred[[k]], type = "l", col = "red",
         xlim = c(0, data.t$tmax),
         ylim = ylimfun(k),
         ylab = ifelse(is.null(ylab), toString(formula(fit$lfit[[k]])[[2]]),
                       ylab[[1]]),
         las = 1,
         xaxt = "n", ...)
    if (k == K) {
      axis(1)
    }
    lines(x = data.t$tk[[k]],
         y = data.t$yk[[k]],
         col = "blue")
    points(x = data.t$tk[[k]],
           y = data.t$yk[[k]],
           pch = 8,
           #xpd = TRUE,
           col = "red")
    if (grid) {
      grid()
    }
    abline(v = data.t$tobs)
  }

  # Axis labels
  mtext(main, 3,
        line = 1, outer = TRUE, font = 2, cex = 1.3)
  mtext(ifelse(is.null(xlab), "Time", xlab), 1,
        line = 2.5, outer = TRUE)

  par(old.par)

}
