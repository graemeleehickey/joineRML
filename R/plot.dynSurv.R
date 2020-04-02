#' Plot a \code{dynSurv} object
#'
#' @description Plots the conditional time-to-event distribution for a
#'   \emph{new} subject calculated using the \code{\link{dynSurv}} function.
#'
#' @param x an object of class \code{dynSurv} calculated by the
#'   \code{\link{dynSurv}} function.
#' @param main an overall title for the plot: see \code{\link[graphics]{title}}.
#' @param xlab a title for the x [time] axis: see \code{\link[graphics]{title}}.
#' @param ylab1 a character vector of the titles for the \emph{K} longitudinal
#'   outcomes y-axes: see \code{\link[graphics]{title}}.
#' @param ylab2 a title for the event-time outcome axis: see
#'   \code{\link[graphics]{title}}.
#' @param grid adds a rectangular grid to an existing plot: see
#'   \code{\link[graphics]{grid}}.
#' @param estimator a character string that can take values \code{mean} or
#'   \code{median} to specify what prediction statistic is plotted from an
#'   objecting inheritting of class \code{dynSurv}. Default is
#'   \code{estimator='median'}. This argument is ignored for non-simulated
#'   \code{dynSurv} objects, i.e. those of \code{type='first-order'}, as in that
#'   case a mode-based prediction is plotted.
#' @param smooth logical: whether to overlay a smooth survival curve (see
#'   \strong{Details}). Defaults to \code{FALSE}.
#' @param ... additional plotting arguments; currently limited to \code{lwd} and
#'   \code{cex}. See \code{\link[graphics]{par}} for details.
#'
#' @details The \code{joineRML} package is based on a semi-parametric model,
#'   such that the baseline hazards function is left unspecified. For
#'   prediction, it might be preferable to have a smooth survival curve. Rather
#'   than changing modelling framework \emph{a prior}, a constrained B-splines
#'   non-parametric median quantile curve is estimated using
#'   \code{\link[cobs]{cobs}}, with a penalty function of \eqn{\lambda=1}, and
#'   subject to constraints of monotonicity and \eqn{S(t)=1}.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords hplot
#' @importFrom cobs cobs
#' @seealso \code{\link{dynSurv}}
#'
#' @references
#'
#' Ng P, Maechler M. A fast and efficient implementation of qualitatively
#' constrained quantile smoothing splines. \emph{Statistical Modelling}. 2007;
#' \strong{7(4)}: 315-328.
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
#' out1 <- dynSurv(fit2, hvd2)
#' plot(out1, main = "Patient 1")
#' }
#'
#' \dontrun{
#' # Monte Carlo simulation with 95% confidence intervals on plot
#'
#' out2 <- dynSurv(fit2, hvd2, type = "simulated", M = 200)
#' plot(out2, main = "Patient 1")
#' }
plot.dynSurv <- function(x, main = NULL, xlab = NULL, ylab1 = NULL,
                         ylab2 = NULL, grid = TRUE, estimator, smooth = FALSE, ...) {

  if (!inherits(x, "dynSurv")) {
    stop("Use only with 'dynSurv' objects.\n")
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
  widths <- widths / sum(widths)
  if (widths[1] < 0.1) {
    widths <- c(0.15, 0.85)
    warning("Longitudinal and event-time time data may be shown on different scales")
  }
  layout(m, widths = widths)
  xticks <- pretty(c(unlist(data.t$tk), pred$u))
  #layout.show(K + 1)
  par(mar = c(0, 4.5, 0, 0), oma = c(4, 0, 3, 0))

  # Fine control plotting arguments
  lwd <- 1
  cex <- 1
  if (!missing(...)) {
    dots <- list(...)
    if ("lwd" %in% names(dots)) {
      lwd <- dots[["lwd"]]
    }
    if ("cex" %in% names(dots)) {
      cex <- dots[["cex"]]
    }
  }

  if (!is.null(ylab1)) {
    if (length(ylab1) != K) {
      stop("Number of longitudinal axes labels does not match number of outcomes.\n")
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
         xaxt = "n",
         ylab = ifelse(is.null(ylab1), toString(formula(fit$lfit[[k]])[[2]]),
                       ylab1[[k]]),
         lwd = lwd)
    points(x = data.t$tk[[k]],
           y = data.t$yk[[k]],
           pch = 8,
           #xpd = TRUE,
           col = "red",
           cex = cex)
    if (k == K) {
      axis(1, at = xticks)
    }
    if (grid) {
      grid()
    }
  }

  ylim <- NULL
  if (x$type == "simulated") {
    ylim <- c(min(pred$lower), max(pred$upper))
  }

  # Plot survival distribution
  par(mar = c(0, 0, 0, 4.5))
  xpts <- pred$u
  if (x$type == "first-order") {
    ypts <- pred$surv
  } else if (x$type == "simulated") {
    if (missing(estimator) || estimator == "median") {
      ypts <- pred$median
    } else if (estimator == "mean") {
      ypts <- pred$mean
    } else {
      stop("estimator must be equal to 'mean' or 'median'\n")
    }
  }
  plot(xpts, ypts,
       yaxt = "n",
       xaxs = "i",
       xlim = c(data.t$tobs, 1.04 * max(pred$u)),
       ylim = ylim,
       type = "s",
       las = 1,
       lwd = lwd)
  if (x$type == "simulated") { # CIs for MC simulated predictions only
    n <- length(xpts)
    xpts2 <- rep(c(xpts, rev(xpts)), each = 2)
    xpts2 <- xpts2[-c(2*n, 2*n + 1)]
    ypts2 <- c(rep(pred$lower, each = 2)[-1], rep(rev(pred$upper), each = 2)[-2*n])
    polygon(x = xpts2, y = ypts2,
            #x = c(xpts, rev(xpts)), y = c(pred$lower, rev(pred$upper)),
            col = "lightgrey",
            border = "lightgrey",
            lwd = 2)
    lines(xpts, ypts, type = "s", lwd = lwd)
  }
  axis(4, las = 1)
  if (grid) {
    grid()
    grid(col = "white")
  }
  abline(v = data.t$tobs, col = "white", lwd = 3, xpd = NA)
  abline(v = data.t$tobs, col = "darkgrey", lty = "dotted", lwd = 3, xpd = FALSE)

  # smoothed survival function
  if (smooth) {
    cobs_fit <- cobs::cobs(xpts, ypts, constraint = "decrease", lambda = 0.1,
                           nknots = 7, pointwise = rbind(c(0, x$data.t$tobs, 1)),
                           print.warn = FALSE, print.mesg = FALSE)
    smooth_pts <- predict(cobs_fit)
    lines(smooth_pts[, 1], smooth_pts[, 2], col = 2, lwd = 2)

    if(x$type == "simulated") {
      # lower CI
      cobs_fit_low <- cobs::cobs(xpts, pred$lower, constraint = "decrease",
                             nknots = 7, pointwise = rbind(c(0, x$data.t$tobs, 1)),
                             print.warn = FALSE, print.mesg = FALSE)
      smooth_pts_low <- predict(cobs_fit_low)
      lines(smooth_pts_low[, 1], smooth_pts_low[, 2], col = 2, lwd = 2, lty = 2)
      # upper CI
      cobs_fit_upp <- cobs::cobs(xpts, pred$upper, constraint = "decrease",
                                 nknots = 7, pointwise = rbind(c(0, x$data.t$tobs, 1)),
                                 print.warn = FALSE, print.mesg = FALSE)
      smooth_pts_upp <- predict(cobs_fit_upp)
      lines(smooth_pts_upp[, 1], smooth_pts_upp[, 2], col = 2, lwd = 2, lty = 2)
    }
  }

  # Axis labels
  mtext(main, 3,
        line = 1, outer = TRUE, font = 2, cex = 1.3)
  mtext(ifelse(is.null(xlab), "Time", xlab), 1,
        line = 2.5, outer = TRUE)
  mtext(ifelse(is.null(ylab2), "Event-free probability", ylab2), 4,
        line = 2.5)

  on.exit(par(old.par))

}
