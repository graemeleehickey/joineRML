#' Dynamic predictions for the time-to-event data sub-model
#'
#' @description Calculates the conditional time-to-event distribution for a
#'   \emph{new} subject from the last observation time given their longitudinal
#'   history data and a fitted \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#' @param newdata a list of \code{data.frame} objects for each longitudinal
#'   outcome for a single new patient in which to interpret the variables named
#'   in the \code{formLongFixed} and \code{formLongRandom} formulae of
#'   \code{object}. As per \code{\link{mjoint}}, the \code{list} structure
#'   enables one to include multiple longitudinal outcomes with different
#'   measurement protocols. If the multiple longitudinal outcomes are measured
#'   at the same time points for each patient, then a \code{data.frame} object
#'   can be given instead of a \code{list}. It is assumed that each data frame
#'   is in long format.
#' @param newSurvData a \code{data.frame} in which to interpret the variables
#'   named in the \code{formSurv} formulae from the \code{mjoint} object. This
#'   is optional, and if omitted, the data will be searched for in
#'   \code{newdata}. Note that no event time or censoring indicator data are
#'   required for dynamic prediction. Defaults to \code{newSurvData=NULL}.
#' @param u an optional time that must be greater than the last observed
#'   measurement time. If omitted (default is \code{u=NULL}), then conditional
#'   failure probabilities are reported for \emph{all} observed failure times in
#'   the \code{mjoint} object data from the last known follow-up time of the
#'   subject.
#' @param horizon an optional horizon time. Instead of specifying a specific
#'   time \code{u} relative to the time origin, one can specify a horizon time
#'   that is relative to the last known follow-up time. The prediction time is
#'   essentially equivalent to \code{horizon} + \eqn{t_{obs}}, where
#'   \eqn{t_{obs}} is the last known follow-up time where the patient had
#'   not yet experienced the event. Default is \code{horizon=NULL}. If
#'   \code{horizon} is non-\code{NULL}, then the output will be reported still
#'   in terms of absolute time (from origin), \code{u}.
#' @param type a character string for whether a first-order
#'   (\code{type="first-order"}) or Monte Carlo simulation approach
#'   (\code{type="simulated"}) should be used for the dynamic prediction.
#'   Defaults to the computationally faster first-order prediction method.
#' @param M for \code{type="simulated"}, the number of simulations to performs.
#'   Default is \code{M=200}.
#' @param scale a numeric scalar that scales the variance parameter of the
#'   proposal distribution for the Metropolis-Hastings algorithm, which
#'   therefore controls the acceptance rate of the sampling algorithm.
#' @param ci a numeric value with value in the interval \eqn{(0, 1)} specifying
#'   the confidence interval level for predictions of \code{type='simulated'}.
#'   If missing, defaults to \code{ci=0.95} for a 95\% confidence interval. If
#'   \code{type='first-order'} is used, then this argument is ignored.
#' @param progress logical: should a progress bar be shown on the console to
#'   indicate the percentage of simulations completed? Default is
#'   \code{progress=TRUE}.
#'
#' @details Dynamic predictions for the time-to-event data sub-model based on an
#'   observed measurement history for the longitudinal outcomes of a new subject
#'   are based on either a first-order approximation or Monte Carlo simulation
#'   approach, both of which are described in Rizopoulos (2011). Namely, given
#'   that the subject was last observed at time \emph{t}, we calculate the
#'   conditional survival probability at time \eqn{u > t} as
#'
#'   \deqn{P[T \ge u | T \ge t; y, \theta] \approx \frac{S(u | \hat{b};
#'   \theta)}{S(t | \hat{b}; \theta)},}
#'
#'   where \eqn{T} is the failure time for the new subject, \eqn{y} is the
#'   stacked-vector of longitudinal measurements up to time \emph{t} and
#'   \eqn{S(u | \hat{b}; \theta)} is the survival function.
#'
#'   \strong{First order predictions}
#'
#'   For \code{type="first-order"}, \eqn{\hat{b}} is the mode
#'   of the posterior distribution of the random effects given by
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
#' @keywords survival
#' @importFrom mvtnorm rmvt
#' @seealso \code{\link{mjoint}}, \code{\link{dynLong}}, and
#'   \code{\link{plot.dynSurv}}.
#'
#' @references
#'
#' Rizopoulos D. Dynamic predictions and prospective accuracy in joint models
#' for longitudinal and time-to-event data. \emph{Biometrics}. 2011;
#' \strong{67}: 819–829.
#'
#' Taylor JMG, Park Y, Ankerst DP, Proust-Lima C, Williams S, Kestin L, et al.
#' Real-time individual predictions of prostate cancer recurrence using joint
#' models. \emph{Biometrics}. 2013; \strong{69}: 206–13.
#'
#' @return A list object inheriting from class \code{dynSurv}. The list returns
#'   the arguments of the function and a \code{data.frame} with first column
#'   (named \code{u}) denoting times and the subsequent columns returning
#'   summary statistics for the conditional failure probabilities For
#'   \code{type="first-order"}, a single column named \code{surv} is appended.
#'   For \code{type="simulated"}, four columns named \code{mean}, \code{median},
#'   \code{lower} and \code{upper} are appended, denoting the mean, median and
#'   lower and upper confidence intervals from the Monte Carlo draws. Additional
#'   objects are returned that are used in the intermediate calculations.
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
#' dynSurv(fit2, hvd2)
#' dynSurv(fit2, hvd2, u = 7) # survival at 7-years only
#'
#' out <- dynSurv(fit2, hvd2, type = "simulated")
#' out
#' }
dynSurv <- function(object, newdata, newSurvData = NULL, u = NULL, horizon = NULL,
                    type = "first-order", M = 200, scale = 2, ci, progress = TRUE) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  if (type != "first-order" && type != "simulated") {
    stop("order must be specified as either 'first-order' or 'simulated")
  }

  data.t <- process_newdata(newdata = newdata,
                            object = object,
                            tobs = NULL,
                            newSurvData = newSurvData)

  if (!is.null(u) && !is.null(horizon)) {
    stop("Cannot specify 'u' and 'horizon' times simultaneously")
  }

  # Construct u (if needed) + check valid
  if (is.null(u)) {
    if (is.null(horizon)) {
      u <- object$dmats$t$tj
      if (length(u) > 1) {
        u <- c(data.t$tobs, u[u > data.t$tobs])
        #u <- u[-length(u)]
      }
    } else {
      u <- horizon + data.t$tobs
    }
  } else { # if u is specified
    if (any(u <= data.t$tobs)) {
      stop("Landmark time(s) must be greater than last observation time")
    }
  }

  # Get posterior mode of [b | data, theta]
  b.hat <- b_mode(data = data.t, theta = object$coefficients)

  if (b.hat$convergence != 0) {
    stop("Optimisation failed")
  }

  #***********************************************
  # First-order prediction
  #***********************************************

  if (type == "first-order") {

    b.hat <- b.hat$par

    # Predicted survival at time t
    S.t <- S(b = b.hat, data = data.t, theta = object$coefficients)

    # Predicted survival at times u
    S.u <- vector(length = length(u))
    for (i in 1:length(u)) {
      data.u <- process_newdata(newdata = newdata,
                                object = object,
                                tobs = u[i],
                                newSurvData = newSurvData)

      S.u[i] <- S(b = b.hat, data = data.u, theta = object$coefficients)
    }

    pred <- data.frame("u" = u, "surv" = S.u / S.t)

  } # end first-order prediction

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
    surv <- matrix(nrow = length(u), ncol = M)
    for (m in 1:M) {
      # Step 1: draw theta
      theta.samp <- thetaDraw(object)
      # Step 2: Metropolis-Hastings simulation
      mh_sim <- b_metropolis(theta.samp, delta.prop, sigma.prop, b.curr, data.t)
      b.curr <- mh_sim$b.curr
      accept <- accept + mh_sim$accept
      # Step 3: predicted survival
      S.t <- S(b = b.curr, data = data.t, theta = theta.samp)
      S.u <- vector(length = length(u))
      for (i in 1:length(u)) {
        data.u <- process_newdata(newdata = newdata,
                                  object = object,
                                  tobs = u[i],
                                  newSurvData = newSurvData)
        S.u[i] <- S(b = b.curr, data = data.u, theta = theta.samp)
      }
      surv[, m] <- S.u / S.t
      if (progress) {
        utils::setTxtProgressBar(pb, m)
      }
    }

    if (progress) {
      close(pb)
    }

    surv.mean <- apply(surv, 1, mean)
    surv.med <- apply(surv, 1, median)
    surv.low <- apply(surv, 1, quantile, prob = alpha)
    surv.upp <- apply(surv, 1, quantile, prob = 1 - alpha)
    pred <- data.frame("u" = u,
                       "mean" = surv.mean,
                       "median" = surv.med,
                       "lower" = surv.low,
                       "upper" = surv.upp)

  } # end simulated prediction

  out <- list("pred" = pred,
              "fit" = object,
              "newdata" = newdata,
              "newSurvData" = newSurvData,
              "u" = u,
              "data.t" = data.t,
              "type" = type,
              "horizon" = horizon,
              "b.hat" = b.hat)

  if (type == "simulated") {
    out$M <- ifelse(type == "simulated", M, NULL)
    out$accept <- ifelse(type == "simulated", accept / M, NULL)
  }

  class(out) <- "dynSurv"
  return(out)

}

