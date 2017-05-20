#' Plot a \code{dynSurv} object
#'
#' @description Plots the conditional time-to-event distribution for a
#'   \emph{new} subject calculated using the \code{\link{dynSurv}} function.
#'
#' @param x an object of class \code{dynSurv} calculated by the
#'   \code{\link{dynSurv}} function.
#' @param main an overall title for the plot: see \code{\link[graphics]{title}}.
#' @param xlab a title for the x [time] axis: see \code{\link[graphics]{title}}.
#' @param ylab1 a title for the \emph{K} longitudinal outcomes y-axes: see
#'   \code{\link[graphics]{title}}.
#' @param ylab2 a title for the event-time outcome axis: see
#'   \code{\link[graphics]{title}}.
#' @param grid adds a rectangular grid to an existing plot: see
#'   \code{\link[graphics]{grid}}.
#' @param ... additional plotting arguments; currently none are used.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords hplot
#' @seealso \code{\link{dynSurv}}
#'
#' @references
#'
#' Rizopoulos D. Dynamic predictions and prospective accuracy in joint models
#' for longitudinal and time-to-event data. \emph{Biometrics}. 2011;
#' \strong{67}: 819–829.
#'
#' @return A plot.
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
#' out <- dynSurv(fit2, hvd2)
#' plot(out, main = "Patient 1")
#' }
plot.dynSurv <- function(x, main = NULL, xlab = NULL, ylab1 = NULL,
                         ylab2 = NULL, grid = TRUE, ...) {

  if (!inherits(x, "dynSurv")) {
    stop("Use only with 'dynSurv' model objects.\n")
  }

  # Extract information we need
  fit <- x$fit
  pred <- x$pred
  newdata <- x$newdata
  newSurvData <- x$newSurvData
  u <- x$u
  data.t <- x$data.t

  # Set-up the plotting grid
  old.par <- par(no.readonly = TRUE)
  K <- fit$dims$K
  m <- cbind(c(1:K), rep(K + 1, K))
  widths <- c(data.t$tobs, max(pred$u) - data.t$tobs)
  layout(m, widths = widths)
  #layout.show(K + 1)
  par(mar = c(0, 4.5, 0, 0), oma = c(4, 0, 3, 0))

  if (!is.null(ylab1)) {
    if (length(ylab1) != K) {
      stop("Number of longitudinal axes labels does not match number of outcomes.")
    }
  }

  # Plot longitudinal markers
  for (k in 1:K) {
    plot(x = data.t$tk[[k]],
         y = data.t$yk[[k]],
         type = "l",
         col = "blue",
         las = 1,
         #xaxs = "i",
         xaxt = ifelse(k == K, "s", "n"),
         ylab = ifelse(is.null(ylab1), toString(formula(fit$lfit[[k]])[[2]]),
                       ylab1[[1]]))
    points(x = data.t$tk[[k]],
           y = data.t$yk[[k]],
           pch = 8,
           #xpd = TRUE,
           col = "red")
    if (grid) {
      grid()
    }

  }

  # Plot survival distribution
  par(mar = c(0, 0, 0, 4.5))
  plot(pred$u, pred$surv,
       yaxt = "n",
       xaxs = "i",
       xlim = c(data.t$tobs, 1.04 * max(pred$u)),
       type = "s",
       las = 1)
  axis(4, las = 1)
  if (grid) {
    grid()
  }
  abline(v = data.t$tobs, col = "white", lwd = 3, xpd = NA)
  abline(v = data.t$tobs, col = "darkgrey", lty = "dotted", lwd = 3, xpd = FALSE)

  # Axis labels
  mtext(main, 3,
        line = 1, outer = TRUE, font = 2, cex = 1.3)
  mtext(ifelse(is.null(xlab), "Time", xlab), 1,
        line = 2.5, outer = TRUE)
  mtext(ifelse(is.null(ylab2), "Event-free probability", ylab), 4,
        line = 2.5)

  par(old.par)

}
