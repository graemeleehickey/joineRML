#' Confidence intervals for model parameters of an \code{mjoint} object
#'
#' @description This function computes confidence intervals for one or more
#'   parameters in a fitted \code{mjoint} object.
#'
#' @param object an object inheriting from class \code{mjoint} for a joint model
#'   of time-to-event and multivariate longitudinal data.
#' @param parm a character string specifying which sub-model parameter
#'   confidence intervals should be returned for. Can be specified as
#'   \code{parm='Longitudinal'} (multivariate longitudinal sub-model),
#'   \code{parm='Event'} (time-to-event sub-model), or \code{parm='both'}
#'   (default).
#' @param level the confidence level required. Default is \code{level=0.95} for
#'   a 95\% confidence interval.
#' @param bootSE an object inheriting from class \code{bootSE} for the
#'   corresponding model. If \code{bootSE=NULL}, the function will attempt to
#'   utilize approximate standard error estimates (if available) calculated from
#'   the empirical information matrix.
#' @param ... additional arguments; currently none are used.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link{mjoint}}, \code{\link{bootSE}}, and
#'   \code{\link[stats]{confint}} for the generic method description.
#'
#' @references
#'
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second
#' Edition. Wiley-Interscience; 2008.
#'
#' Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
#' measurements and event time data. \emph{Biostatistics.} 2000; \strong{1(4)}:
#' 465-480.
#'
#' Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the joint
#' analysis of time-to-event and multiple longitudinal variables. \emph{Stat
#' Med.} 2002; \strong{21}: 2369-2382.
#'
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
#' measured with error. \emph{Biometrics.} 1997; \strong{53(1)}: 330-339.
#'
#' @return A matrix containing the confidence intervals for either the
#'   longitudinal, time-to-event, or both sub-models.
#' @export
#'
#' @examples
#' # Fit a classical univariate joint model with a single longitudinal outcome
#' # and a single time-to-event outcome
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#'
#' gamma <- c(0.1059417, 1.0843359)
#' sigma2 <- 0.03725999
#' beta <- c(4.9988669999, -0.0093527634, 0.0004317697)
#' D <- matrix(c(0.128219108, -0.006665505, -0.006665505, 0.002468688),
#'             nrow = 2, byrow = TRUE)
#'
#' set.seed(1)
#' fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
#'     formLongRandom = ~ time | num,
#'     formSurv = Surv(fuyrs, status) ~ age,
#'     data = hvd,
#'     timeVar = "time",
#'     inits = list(gamma = gamma, sigma2 = sigma2, beta = beta, D = D),
#'     control = list(nMCscale = 2, burnin = 5)) # controls for illustration only
#'
#' confint(fit1, parm = "Longitudinal")
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
#'     verbose = TRUE)
#' confint(fit2)
#' }
confint.mjoint <- function(object, parm = c("Both", "Longitudinal", "Event"),
                           level = 0.95, bootSE = NULL, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  parm <- match.arg(parm)

  dims <- object$dims
  num.d <- with(dims, sum(r) * (sum(r) + 1) / 2) # number of RE covariance params
  num.b <- sum(dims$p)                           # number of beta coefficients
  num.s <- dims$K                                # number of error variances
  num.g <- with(dims, q + K)                     # number of gamma coefficients

  q <- qnorm((1 - level) / 2, lower.tail = FALSE)

  beta <- object$coefficients$beta
  beta.inds <- (num.d + 1):(num.d + num.b)
  if (is.null(bootSE)) {
    beta.se <- sqrt(diag(vcov(object)))[beta.inds]
  } else {
    beta.se <- bootSE$beta.se
  }

  beta.ci <- cbind("lower" = beta - q * beta.se,
                   "upper" = beta + q * beta.se)

  gamma <- object$coefficients$gamma
  gamma.inds <- (num.d + num.b + num.s + 1):(num.d + num.b + num.s + num.g)
  if (is.null(bootSE)) {
    gamma.se <- sqrt(diag(vcov(object)))[gamma.inds]
  } else {
    gamma.se <- bootSE$gamma.se
  }

  gamma.ci <- cbind("lower" = gamma - q * gamma.se,
                    "upper" = gamma + q * gamma.se)

  if (parm == "Both") {
    ci.mat <- rbind(beta.ci, gamma.ci)
  } else if (parm == "Longitudinal") {
    ci.mat <- beta.ci
  } else {
    ci.mat <- gamma.ci
  }

  colnames(ci.mat) <- c(paste0((1 - level)*50, "%"),
                        paste0(100 - (1 - level)*50, "%"))
  return(ci.mat)

}
