#' Dynamic predictions for the longitudinal data sub-model
#'
#' @description Calculates the conditional expected longitudinal values for a
#'   \emph{new} subject from the last observation time given their longitudinal
#'   history data and a fitted \code{mjoint} object.
#'
#' @inheritParams dynSurv
#' @param ntimes an integer controlling the number of points to discretize the
#'   extrapolated time region into. Default is \code{ntimes=100}.
#' @param level an optional integer giving the level of grouping to be used in
#'   extracting the residuals from object. Level values increase from outermost
#'   to innermost grouping, with level 0 corresponding to the population model
#'   fit and level 1 corresponding to subject-specific model fit. Defaults to
#'   \code{level=1}.
#'
#' @details Dynamic predictions for the longitudinal data sub-model based on an
#'   observed measurement history for the longitudinal outcomes of a new subject
#'   are based on either a first-order approximation or Monte Carlo simulation
#'   approach, both of which are described in Rizopoulos (2011). Namely, given
#'   that the subject was last observed at time \emph{t}, we calculate the
#'   conditional expectation of each longitudinal outcome at time \emph{u} as
#'
#'   \deqn{E[y_k(u) | T \ge t, y, \theta] \approx x^T(u)\beta_k +
#'   z^T(u)\hat{b}_k,}
#'
#'   where \eqn{T} is the failure time for the new subject, and \eqn{y} is the
#'   stacked-vector of longitudinal measurements up to time \emph{t}.
#'
#'   \strong{First order predictions}
#'
#'   For \code{type="first-order"}, \eqn{\hat{b}} is the mode of the posterior
#'   distribution of the random effects given by
#'
#'   \deqn{\hat{b} = {\arg \max}_b f(b | y, T \ge t; \theta).}
#'
#'   The predictions are based on plugging in \eqn{\theta = \hat{\theta}}, which
#'   is extracted from the \code{mjoint} object.
#'
#'   \strong{Monte Carlo simulation predictions}
#'
#'   For \code{type="simulated"}, \eqn{\theta} is drawn from a multivariate
#'   normal distribution with means \eqn{\hat{\theta}} and variance-covariance
#'   matrix both extracted from the fitted \code{mjoint} object via the
#'   \code{coef()} and \code{vcov()} functions. \eqn{\hat{b}} is drawn from the
#'   the posterior distribution of the random effects
#'
#'   \deqn{f(b | y, T \ge t; \theta)}
#'
#'   by means of a Metropolis-Hasting algorithm with independent multivariate
#'   non-central \emph{t}-distribution proposal distributions with
#'   non-centrality parameter \eqn{\hat{b}} from the first-order prediction and
#'   variance-covariance matrix equal to \code{scale} \eqn{\times} the inverse
#'   of the negative Hessian of the posterior distribution. The choice os
#'   \code{scale} can be used to tune the acceptance rate of the
#'   Metropolis-Hastings sampler. This simulation algorithm is iterated \code{M}
#'   times, at each time calculating the conditional survival probability.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords multivariate
#' @seealso \code{\link{mjoint}}, \code{\link{dynSurv}}.
#'
#' @references
#'
#' Rizopoulos D. Dynamic predictions and prospective accuracy in joint models
#' for longitudinal and time-to-event data. \emph{Biometrics}. 2011;
#' \strong{67}: 819–829.
#'
#' @return A list object inheriting from class \code{dynLong}. The list returns
#'   the arguments of the function and a list containing \emph{K}
#'   \code{data.frame}s of 2 columns, with first column (named
#'   \code{timeVar[k]}; see \code{\link{mjoint}}) denoting times and the second
#'   column (named \code{y.pred}) denoting the expected outcome at each time
#'   point.
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
#' dynLong(fit2, hvd2)
#' dynLong(fit2, hvd2, u = 7) # outcomes at 7-years only
#'
#' out <- dynLong(fit2, hvd2, type = "simulated")
#' out
#' }
dynLong <- function(object, newdata, newSurvData = NULL, u = NULL,
                    type = "first-order", M = 200, scale = 1.6, ci,
                    progress = TRUE, ntimes = 100, level = 1) {

  K <- object$dims$K
  p <- object$dims$p
  r <- object$dims$r
  beta.inds <- cumsum(c(0, p))
  b.inds <- cumsum(c(0, r))
  beta <- object$coefficients$beta

  if (class(newdata) != "list") {
    balanced <- TRUE
    newdata <- list(newdata)
    if (K > 1) {
      for (k in 2:K) {
        newdata[[k]] <- newdata[[1]]
      }
    }
  } else {
    balanced <- (length(unique(newdata)) == 1)
  }
  if (length(newdata) != K) {
    stop(paste("The number of datasets expected is K =", K))
  }

  data.t <- process_newdata(newdata = newdata,
                            object = object,
                            tobs = NULL,
                            newSurvData = newSurvData)

  if (is.null(u)) {
    pred.times <- seq(from = data.t$tobs, to = data.t$tmax,
                      length.out = ntimes)
  } else{
    if (length(u) > 1) {
      stop("Only a single time can be accepted.\n")
    }
    if (u < data.t$tobs) {
      stop("Predictions can only be made after the last observation time.\n")
    }
    pred.times <- u
    ntimes <- 1
  }

  # Get posterior mode of [b | data, theta]
  b.hat <- b_mode(data = data.t, theta = object$coefficients)

  if (b.hat$convergence != 0) {
    stop("Optimisation failed")
  }

  # Extrapolated data
  newdata2 <- list()
  for (k in 1:K) {
    timecol <- which(names(newdata[[k]]) == object$timeVar[[k]])
    newdata2[[k]] <- newdata[[k]][, -timecol]
    baseline <- newdata2[[k]][nrow(newdata2[[k]]), , drop = FALSE][rep(1, ntimes), ]
    newdata2[[k]] <- cbind(baseline, pred.times)
    names(newdata2[[k]])[ncol(newdata2[[k]])] <- object$timeVar[[k]]
  }

  # Design matrices for extrapolated data
  # NB: ignore everything else in data.t2
  data.t2 <- process_newdata(newdata = newdata2,
                             object = object,
                             tobs = NULL,
                             newSurvData = newSurvData)

  #***********************************************
  # First-order prediction
  #***********************************************

  if (type == "first-order") {

    # Predicted y over extrapolated data
    y.pred <- list()
    for (k in 1:K) {
      beta.k <- beta[(beta.inds[k] + 1):(beta.inds[k + 1])]
      Xbeta.k <- data.t2$Xk[[k]] %*% beta.k
      y.pred[[k]] <- Xbeta.k
      if (level == 1) {
        Zb.k <- data.t2$Zk[[k]] %*% b.hat$par[(b.inds[k] + 1):(b.inds[k + 1])]
        y.pred[[k]] <- y.pred[[k]] + Zb.k
      }
    }

    pred <- lapply(y.pred, function(x) {
      df <- cbind(pred.times, x)
      df <- as.data.frame(df)
      names(df) <- c("time", "y.pred")
      df
    })
  }

  #***********************************************
  # Simulated prediction
  #***********************************************

  if (type == "simulated") {

    if (progress) {
      cat("\n")
      pb <- utils::txtProgressBar(min = 0, max = M, style = 3)
    }

    if (missing(ci)) {
      ci <- 0.95
    } else {
      if (!is.numeric(ci) || (ci < 0) || (ci > 1)) {
        stop("ci argument has been misspecified.\n")
      }
    }
    alpha <- (1 - ci) / 2

    b.curr <- delta.prop <- b.hat$par
    sigma.prop <- -solve(b.hat$hessian) * scale
    accept <- 0
    y.pred <- list()
    for (k in 1:K) {
      y.pred[[k]] <- matrix(nrow = ntimes, ncol = M)
    }
    for (m in 1:M) {
      # Step 1: draw theta
      theta.samp <- thetaDraw(object)
      # Step 2: Metropolis-Hastings simulation
      mh_sim <- b_metropolis(theta.samp, delta.prop, sigma.prop, b.curr, data.t)
      b.curr <- mh_sim$b.curr
      accept <- accept + mh_sim$accept
      # Step 3: predicted longitudinal value
      for (k in 1:K) {
        beta.k <- theta.samp$beta[(beta.inds[k] + 1):(beta.inds[k + 1])]
        Xbeta.k <- data.t2$Xk[[k]] %*% beta.k
        y.pred[[k]][, m] <- Xbeta.k
        if (level == 1) {
          Zb.k <- data.t2$Zk[[k]] %*% b.curr[(b.inds[k] + 1):(b.inds[k + 1])]
          y.pred[[k]][, m] <- y.pred[[k]][, m] + Zb.k
        }
      }
      if (progress) {
        utils::setTxtProgressBar(pb, m)
      }
    }

    if (progress) {
      close(pb)
    }

    pred <- lapply(y.pred, function(x) {
      long.mean <- apply(x, 1, mean)
      long.med <- apply(x, 1, median)
      long.se <- apply(x, 1, sd)
      long.low <- apply(x, 1, quantile, prob = alpha)
      long.upp <- apply(x, 1, quantile, prob = 1 - alpha)
      df <- cbind(pred.times, long.mean, long.med, long.se, long.low, long.upp)
      df <- as.data.frame(df)
      names(df) <- c("time", "mean", "median", "se", "lower", "upper")
      df
    })

  }

  names(pred) <- names(object$formLongFixed)
  out <- list("pred" = pred,
              "fit" = object,
              "newdata" = newdata,
              "newSurvData" = newSurvData,
              "u" = u,
              "data.t" = data.t,
              "type" = type)

  if (type == "simulated") {
    out$M <- ifelse(type == "simulated", M, NULL)
    out$accept <- ifelse(type == "simulated", accept / M, NULL)
  }

  class(out) <- "dynLong"
  return(out)

}
