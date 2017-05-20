#' Dynamic predictions for the time-to-event submodel
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
#'   failure probabilities are reported for all observed failure times in the
#'   \code{mjoint} object data.
#'
#' @details Dynamic predictions for the time-to-event submodel based on an
#'   observed measurement history for the longitudinal outcomes of a new subject
#'   are based on a first-order approximation described in Rizopoulos (2011).
#'   Namely, given that the subject was last observed at time \emph{t}, we
#'   calculate the conditional survival probability at time \eqn{u > t} as
#'
#'   \deqn{P[T \geq u | T \geq u; y, \theta] \approx
#'   \frac{S(u | \hat{b}; \theta)}{S(t | \hat{b}; \theta)},}
#'
#'   where \eqn{T} is the failure time for the new subject, \eqn{y} is the
#'   stacked-vector of longitudinal measurements, \eqn{S(u | \hat{b}; \theta)}
#'   is the survival function, and \eqn{\hat{b}} is the mode of the posterior
#'   distribution of the random effects given by
#'
#'   \deqn{\hat{b} = {\arg \max}_b f(b | y, T \geq t; \theta).}
#'
#'   The predictions are based on plugging in \eqn{\theta = \hat{\theta}}, which
#'   is extracted from the \code{mjoint} object.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords survival
#' @seealso \code{\link{mjoint}}
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
#'   the arguments of the function and a \code{data.frame} of 2 columns, with
#'   first column (named \code{u}) denoting times and the second column (named
#'   \code{surv}) denoting the conditional failure probability.
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
#' }
dynSurv <- function(object, newdata, newSurvData = NULL, u = NULL) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  data.t <- process_newdata(newdata = newdata,
                            object = object,
                            tobs = NULL,
                            newSurvData = newSurvData)

  # Construct u (if needed) + check valid
  if (is.null(u)) {
    u <- object$dmats$t$tj
    if (length(u) > 1) {
      u <- c(data.t$tobs, u[u > data.t$tobs])
      #u <- u[-length(u)]
    }
  } else {
    if (any(u <= data.t$tobs)) {
      stop("Landmark time(s) must be greater than last observation time")
    }
  }

  # Get posterior mode of [b | data, theta]
  b.hat <- b.mode(data = data.t, theta = object$coefficients)

  if (b.hat$convergence == 0) {
    b.hat <- b.hat$par
  } else {
    stop("Optimisation failed")
  }

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

  out <- list("pred" = pred,
              "fit" = object,
              "newdata" = newdata,
              "newSurvData" = newSurvData,
              "u" = u,
              "data.t" = data.t)

  class(out) <- "dynSurv"
  return(out)

}

