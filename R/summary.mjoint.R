#' Summary of an \code{mjoint} object
#'
#' @description This function provides a summary of an \code{mjoint} object.
#'
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link{mjoint}}, \code{\link{mjoint.object}}, and
#'   \code{\link[base]{summary}} for the generic method description.
#'
#' @references
#'
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
#' measured with error. \emph{Biometrics.} 1997; \strong{53(1)}: 330-339.
#'
#' Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
#' measurements and event time data. \emph{Biostatistics.} 2000; \strong{1(4)}:
#' 465-480.
#'
#' Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the joint
#' analysis of time-to-event and multiple longitudinal variables. \emph{Stat
#' Med.} 2002; \strong{21}: 2369-2382.
#'
#' @return A list containing the coefficient matrices for the longitudinal and
#'   time-to-event sub-models; variance-covariance matrix for the random
#'   effects; residual error variances; log-likelihood of joint model; AIC and
#'   BIC statistics; and model fit objects.
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
#' summary(fit2)
#' }
summary.mjoint <- function(object, bootSE = NULL, ...) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  dims <- object$dims
  num.d <- with(dims, sum(r) * (sum(r) + 1) / 2) # number of RE covariance params
  num.b <- sum(dims$p)                           # number of beta coefficients
  num.s <- dims$K                                # number of error variances
  num.g <- with(dims, q + K)                     # number of gamma coefficients

  beta <- object$coefficients$beta
  beta.inds <- (num.d + 1):(num.d + num.b)
  if (is.null(bootSE) & !is.null(object$Hessian)) {
    beta.se <- sqrt(diag(vcov(object)))[beta.inds]
  } else if (!is.null(bootSE)) {
    beta.se <- bootSE$beta.se
  } else {
    beta.se <- rep(NA, length(beta))
  }
  coefs.beta <- cbind(
    "Value" = beta,
    "Std.Err" = beta.se,
    "z-value" = beta / beta.se,
    "p-value" = 2 * pnorm(abs(beta / beta.se), lower.tail = FALSE))

  gamma <- object$coefficients$gamma
  gamma.inds <- (num.d + num.b + num.s + 1):(num.d + num.b + num.s + num.g)
  if (is.null(bootSE) & !is.null(object$Hessian)) {
    gamma.se <- sqrt(diag(vcov(object)))[gamma.inds]
  } else if (!is.null(bootSE)) {
    gamma.se <- bootSE$gamma.se
  } else {
    gamma.se <- rep(NA, length(gamma))
  }
  coefs.gamma <- cbind(
    "Value" = gamma,
    "Std.Err" = gamma.se,
    "z-value" = gamma / gamma.se,
    "p-value" = 2 * pnorm(abs(gamma / gamma.se), lower.tail = FALSE))

  out <- list("coefs.long" = coefs.beta,
              "coefs.surv" = coefs.gamma,
              D = getVarCov(object),
              sigma = sqrt(object$coefficients$sigma2),
              logLik = as.vector(logLik(object)),
              AIC = AIC(object),
              BIC = AIC(object, k = log(dims$n)))

  out$formLongFixed <- object$formLongFixed
  out$formLongRandom <- object$formLongRandom
  out$formSurv <- object$formSurv
  out$sfit <- object$sfit
  out$lfit <- object$lfit
  out$timeVar <- object$timeVar
  out$dims <- object$dims
  out$control <- object$control
  out$finalnMC <- object$finalnMC
  out$call <- object$call
  out$comp.time <- object$comp.time
  out$conv <- object$conv
  out$se.type <- ifelse(is.null(bootSE) & !is.null(object$Hessian), "approx",
                        ifelse(!is.null(bootSE), "boot", "none"))
  if (!is.null(bootSE)) {
    out$boot.time <- bootSE$boot.time
    out$nboot <- bootSE$nboot
  }

  class(out) <- "summary.mjoint"
  out

}
