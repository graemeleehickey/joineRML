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
#'   approach. See Rizopoulos (2011) for the prediction framework. Given
#'   that the subject was last observed at time \emph{t}, we calculate the
#'   conditional expectation of each longitudinal outcome at time \emph{u} as
#'
#'   \deqn{E[y_k(u) | T \ge t, y, \theta] \approx x^T(u)\beta_k +
#'   z^T(u)\hat{b}_k,}
#'
#'   where \eqn{T} is the failure time for the new subject, and \eqn{y} is the
#'   stacked-vector of longitudinal measurements up to time \emph{t}.
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
#'   For simulated predictions, \code{proposals} records the number of proposals
#'   for each of the \code{M} accepted random-effects draws, \code{accept} is
#'   \code{M / sum(proposals)}, and \code{sampler} is
#'   \code{"conditional-rejection"}. The former \code{scale} setting is ignored.
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
                    progress = TRUE, ntimes = 100, level = 1,
                    max.proposals = 10000) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.")
  }
  if (length(type) != 1L || is.na(type) || !type %in% c("first-order", "simulated")) {
    stop("type must be specified as either 'first-order' or 'simulated'")
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
  K <- object$dims$K
  p <- object$dims$p
  r <- object$dims$r
  beta.inds <- cumsum(c(0, p))
  b.inds <- cumsum(c(0, r))
  beta <- object$coefficients$beta

  if (!isa(newdata, "list")) {
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

  # The conditional mode is needed only for first-order predictions.
  b.hat <- NULL
  if (type == "first-order") {
    b.hat <- b_mode(data = data.t, theta = object$coefficients)
    if (b.hat$convergence != 0) stop("Optimisation failed")
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
      on.exit(close(pb), add = TRUE)
    }

    alpha <- (1 - ci) / 2

    proposals <- integer(M)
    y.pred <- list()
    for (k in 1:K) {
      y.pred[[k]] <- matrix(nrow = ntimes, ncol = M)
    }
    for (m in seq_len(M)) {
      # Draw parameters, then random effects conditional on those parameters.
      draw <- .mjoint_prediction_draw(object, data.t, max.proposals, m, M)
      theta.samp <- draw$theta
      b.curr <- draw$b
      proposals[m] <- draw$proposals
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
    out$M <- M
    out$sampler <- "conditional-rejection"
    out$proposals <- proposals
    out$max.proposals <- max.proposals
    out$accept <- M / sum(proposals)
  }

  class(out) <- "dynLong"
  return(out)

}
