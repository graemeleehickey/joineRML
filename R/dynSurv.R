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
#'   survival probabilities are reported for \emph{all} observed failure times in
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
#'   Supplied horizons must be numeric, finite, and strictly positive.
#' @param type a character string for whether a first-order
#'   (\code{type="first-order"}) or Monte Carlo simulation approach
#'   (\code{type="simulated"}) should be used for the dynamic prediction.
#'   Defaults to the computationally faster first-order prediction method.
#' @param M for \code{type="simulated"}, the number of simulations. Must be a positive integer.
#'   Default is \code{M=200}.
#' @param scale deprecated; ignored. Supplying it produces a warning. The
#'   conditional rejection sampler does not require a proposal scale.
#' @param max.proposals maximum number of proposals per parameter draw for
#'   \code{type="simulated"}. Must be a positive integer; defaults to
#'   \code{10000}. The prediction stops if this limit is exhausted.
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
#'   approach. See Rizopoulos (2011) for the prediction framework. Given
#'   that the subject was last observed at time \emph{t}, we calculate the
#'   conditional survival probability at time \eqn{u > t} as
#'
#'   \deqn{P[T \ge u | T \ge t; y, \theta] \approx \frac{S(u | \hat{b};
#'   \theta)}{S(t | \hat{b}; \theta)},}
#'
#'   where \eqn{T} is the failure time for the new subject, \eqn{y} is the
#'   stacked-vector of longitudinal measurements up to time \emph{t} and
#'   \eqn{S(u | \hat{b}; \theta)} is the survival function. The ratio is
#'   evaluated from hazard increments between \eqn{t} and \eqn{u}, avoiding
#'   division of survival probabilities that may both round to zero.
#'
#'   Prediction matrices use the transformations, factor levels, and contrasts
#'   retained when the model was fitted. Changing \code{options("contrasts")}
#'   does not change predictions from that fit. New factor levels are rejected.
#'   Older fits are checked against their stored random-effects matrices;
#'   inconsistent bases at measurement and event times require refitting with
#'   \code{mjoint()}. Transformations that depend on the data must support R's
#'   \code{\link[stats]{makepredictcall}} mechanism or use fixed constants.
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
#'   For \code{type="simulated"}, the non-hazard parameters are drawn jointly
#'   from a multivariate normal approximation with residual variances on the
#'   log scale. Write \eqn{v_k = \sigma_k^2} and let \eqn{g(\theta)} replace
#'   each \eqn{v_k} by \eqn{\log(v_k)}, leaving the other parameters unchanged.
#'   The mean is \eqn{g(\hat{\theta})}, and the covariance matrix is
#'   \eqn{J V J^T}, where \eqn{V} is returned by \code{vcov(object)} and
#'   \eqn{J} is diagonal, with entries \eqn{1 / \hat{v}_k} for the residual
#'   variances and 1 elsewhere. This delta-method transformation includes the
#'   covariances with the other parameters. Exponentiating the log-variance
#'   draws gives positive residual variances, with medians equal to their fitted
#'   estimates. Their original-scale means and covariances are not preserved
#'   exactly. Fitted estimates and reported standard errors retain their
#'   original scales. Invalid fitted variances, unusable covariance matrices,
#'   and numerical overflow or underflow in the draws produce an error.
#'
#'   This is an asymptotic approximation; positivity alone does not establish
#'   its accuracy in small samples or near a variance boundary. Random-effects
#'   covariance matrices and baseline-hazard increments retain their existing
#'   adjustments to ensure positive definiteness and positivity, respectively.
#'
#'   For each parameter draw, random effects are sampled from
#'   \eqn{f(b | y, T \ge t; \theta)} by Gaussian rejection sampling. The
#'   longitudinal conditional distribution is \eqn{N(m,V)}. Its density is
#'   multiplied by the survival contribution \eqn{\exp\{\ell(b)\}}, where
#'   \eqn{\ell(b) = -\sum_j h_j \exp(a_j^T b)} is concave and \eqn{h_j \ge 0}.
#'   At an approximate conditional mode \eqn{b_0}, put
#'   \eqn{g_0 = \nabla\ell(b_0)}. A tangent bound gives the proposal
#'   \eqn{N(m + V g_0,V)}, accepted with log probability
#'   \eqn{\ell(b) - \ell(b_0) - g_0^T(b-b_0)}. This samples the conditional
#'   distribution without burn-in; an approximate mode affects efficiency,
#'   rather than the target distribution. Acceptance is evaluated on the log
#'   scale to avoid cancellation and overflow.
#'
#'   Each parameter draw is held fixed until a random-effects draw is accepted.
#'   Exhausting \code{max.proposals} stops the prediction and identifies the
#'   simulation and conditioning time; parameter draws are never discarded and
#'   replaced. Strong survival conditioning, particularly in several independent
#'   random-effects directions, can require many proposals. Increasing
#'   \code{max.proposals} allows more attempts but does not improve the
#'   approximation used to sample the model parameters. The conditional sampler
#'   alone does not establish prediction-interval coverage.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords survival
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
#'   summary statistics for the conditional survival probabilities. For
#'   \code{type="first-order"}, a single column named \code{surv} is appended.
#'   For \code{type="simulated"}, four columns named \code{mean}, \code{median},
#'   \code{lower} and \code{upper} are appended, denoting the mean, median and
#'   lower and upper confidence limits from the Monte Carlo draws. Additional
#'   objects are returned that are used in the intermediate calculations.
#'   For simulated predictions, \code{proposals} records the number of proposals
#'   for each of the \code{M} accepted random-effects draws, \code{accept} is
#'   \code{M / sum(proposals)}, and \code{sampler} is
#'   \code{"conditional-rejection"}. The former \code{scale} setting is ignored.
#'   The \code{b.hat} component contains the conditional mode for first-order
#'   predictions and is \code{NULL} for simulated predictions.
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
                    type = "first-order", M = 200, scale = 2, ci, progress = TRUE,
                    max.proposals = 10000) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  if (length(type) != 1L || is.na(type) || !type %in% c("first-order", "simulated")) {
    stop("type must be specified as either 'first-order' or 'simulated")
  }

  if (!missing(scale)) {
    warning("'scale' is deprecated and ignored; simulated predictions use conditional rejection sampling.",
            call. = FALSE)
  }
  if (type == "simulated") {
    if (missing(ci)) ci <- 0.95
    .mjoint_simulation_controls(M, max.proposals, ci)
  }

  object$designs <- .mjoint_prediction_designs(object)
  data.t <- process_newdata(newdata = newdata,
                            object = object,
                            tobs = NULL,
                            newSurvData = newSurvData)

  if (!is.null(u) && !is.null(horizon)) {
    stop("Cannot specify 'u' and 'horizon' times simultaneously")
  }

  if (!is.null(horizon) &&
      (!is.numeric(horizon) || length(horizon) == 0L ||
       any(!is.finite(horizon)) || any(horizon <= 0))) {
    stop("'horizon' must contain finite, strictly positive numeric values")
  }

  # Construct u (if needed) + check valid
  if (is.null(u)) {
    if (is.null(horizon)) {
      u <- object$dmats$t$tj
      u <- c(data.t$tobs, u[u > data.t$tobs])
    } else {
      u <- horizon + data.t$tobs
    }
  } else { # if u is specified
    if (!is.numeric(u) || length(u) == 0L || any(!is.finite(u))) {
      stop("'u' must contain finite numeric values")
    }
    if (any(u <= data.t$tobs)) {
      stop("Landmark time(s) must be greater than last observation time")
    }
  }

  # Event-time designs are shared by all parameter draws.
  data.u <- lapply(u, function(time) {
    process_newdata(newdata = newdata, object = object, tobs = time,
                    newSurvData = newSurvData)
  })

  # The conditional mode is needed only for first-order predictions.
  b.hat <- NULL
  if (type == "first-order") {
    b.hat <- b_mode(data = data.t, theta = object$coefficients)
    if (b.hat$convergence != 0) stop("Optimisation failed")
  }

  #***********************************************
  # First-order prediction
  #***********************************************

  if (type == "first-order") {

    b.hat <- b.hat$par

    survival <- vapply(data.u, function(d) {
      S(b = b.hat, data = d, theta = object$coefficients, start = data.t$tj.ind)
    }, numeric(1))
    pred <- data.frame("u" = u, "surv" = survival)

  } # end first-order prediction

  #***********************************************
  # Simulated prediction
  #***********************************************

  if (type == "simulated") {

    if (progress) {
      cat("\n")
      pb <- utils::txtProgressBar(min = 0, max = M, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    alpha <- (1 - ci) / 2

    proposals <- integer(M)
    surv <- matrix(nrow = length(u), ncol = M)
    for (m in seq_len(M)) {
      # Draw parameters, then random effects conditional on those parameters.
      draw <- .mjoint_prediction_draw(object, data.t, max.proposals, m, M)
      theta.samp <- draw$theta
      b.curr <- draw$b
      proposals[m] <- draw$proposals
      # Accumulate hazards only after the conditioning time.
      surv[, m] <- vapply(data.u, function(d) {
        S(b = b.curr, data = d, theta = theta.samp, start = data.t$tj.ind)
      }, numeric(1))
      if (progress) {
        utils::setTxtProgressBar(pb, m)
      }
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
    out$M <- M
    out$sampler <- "conditional-rejection"
    out$proposals <- proposals
    out$max.proposals <- max.proposals
    out$accept <- M / sum(proposals)
  }

  class(out) <- "dynSurv"
  return(out)

}

