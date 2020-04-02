#' Simulate data from a joint model
#'
#' @description This function simulates multivariate longitudinal and
#'   time-to-event data from a joint model.
#'
#' @param n the number of subjects to simulate data for.
#' @param ntms the maximum number of (discrete) time points to simulate repeated
#'   longitudinal measurements at.
#' @param beta a matrix of \code{dim=c(K,4)} specifying the coefficients of the
#'   fixed effects. The order in each row is intercept, time, a continuous
#'   covariate, and a binary covariate.
#' @param gamma.x a vector of \code{length=2} specifying the coefficients for
#'   the time-to-event baseline covariates, in the order of a continuous
#'   covariate and a binary covariate.
#' @param gamma.y a vector of \code{length=K} specifying the latent association
#'   parameters for each longitudinal outcome.
#' @param sigma2 a vector of \code{length=K} specifying the residual standard
#'   errors.
#' @param D a positive-definite matrix specifying the variance-covariance
#'   matrix. If \code{model='int'}, the matrix has dimension \code{dim=c(K, K)},
#'   else if \code{model='intslope'}, the matrix has dimension \code{dim =c(2K,
#'   2K)}. If \code{D=NULL} (default), an identity matrix is assumed.
#' @param df a non-negative scalar specifying the degrees of freedom for the
#'   random effects if sampled from a multivariate \emph{t}-distribution. The
#'   default is \code{df=Inf}, which corresponds to a multivariate normal
#'   distribution.
#' @param model follows the model definition in the \code{\link[joineR]{joint}}
#'   function. See \strong{Details} for choices.
#' @param theta0,theta1 parameters controlling the failure rate. See Details.
#' @param censoring logical: if \code{TRUE}, includes an independent censoring
#'   time.
#' @param censlam a scale (\eqn{> 0}) parameter for an exponential distribution
#'   used to simulate random censoring times for when \code{censoring=TRUE}.
#' @param truncation logical: if \code{TRUE}, adds a truncation time for a
#'   maximum event time.
#' @param trunctime a truncation time for use when \code{truncation=TRUE}.
#'
#' @details The function \code{simData} simulates data from a joint model,
#'   similar to that performed in Henderson et al. (2000). It works by first
#'   simulating multivariate longitudinal data for all possible follow-up times
#'   using random draws for the multivariate Gaussian random effects and
#'   residual error terms. Data can be simulated assuming either
#'   random-intercepts only in each of the longitudinal sub-models, or
#'   random-intercepts and random-slopes. Currently, all models must have the
#'   same structure. The failure times are simulated from proportional hazards
#'   time-to-event models using the following methodologies:
#'
#'   \describe{
#'
#'   \item{\code{model="int"}}{The baseline hazard function is specified to be
#'   an exponential distribution with
#'
#'   \deqn{\lambda_0(t) = \exp{\theta_0}.}
#'
#'   Simulation is conditional on known time-independent effects, and the
#'   methodology of Bender et al. (2005) is used to simulate the failure time.}
#'
#'   \item{\code{model="intslope"}}{The baseline hazard function is specified to
#'   be a Gompertz distribution with
#'
#'   \deqn{\lambda_0(t) = \exp{\theta_0 + \theta_1 t}.}
#'
#'   In the usual representation of the Gompertz distribution, \eqn{\theta_1} is
#'   the shape parameter, and the scale parameter is equivalent to
#'   \eqn{\exp(\theta_0)}. Simulation is conditional on on a predictable
#'   (linear) time-varying process, and the methodology of Austin (2012) is used
#'   to simulate the failure time.}
#'
#'   }
#'
#' @author Pete Philipson (\email{peter.philipson1@@newcastle.ac.uk}) and Graeme
#'   L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords datagen multivariate survival
#'
#' @references
#'
#' Austin PC. Generating survival times to simulate Cox proportional hazards
#' models with time-varying covariates. \emph{Stat Med.} 2012; \strong{31(29)}:
#' 3946-3958.
#'
#' Bender R, Augustin T, Blettner M. Generating survival times to simulate Cox
#' proportional hazards models. \emph{Stat Med.} 2005; \strong{24}: 1713-1723.
#'
#' Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
#' measurements and event time data. \emph{Biostatistics.} 2000; \strong{1(4)}:
#' 465-480.
#'
#' @importFrom mvtnorm rmvt
#'
#' @return A list of 2 \code{data.frame}s: one recording the requisite
#'   longitudinal outcomes data, and one recording the time-to-event data.
#' @export
#'
#' @examples
#' beta <- rbind(c(0.5, 2, 1, 1),
#' c(2, 2, -0.5, -1))
#' D <- diag(4)
#' D[1, 1] <- D[3, 3] <- 0.5
#' D[1, 2] <- D[2, 1] <- D[3, 4] <- D[4, 3] <- 0.1
#' D[1, 3] <- D[3, 1] <- 0.01
#'
#' sim <- simData(n = 250, beta = beta, D = D, sigma2 = c(0.25, 0.25),
#'                censlam = exp(-0.2), gamma.y = c(-.2, 1), ntms = 8)
simData <- function(n = 100, ntms = 5, beta = rbind(c(1, 1, 1, 1), c(1, 1, 1, 1)),
                    gamma.x = c(1, 1), gamma.y = c(0.5, -1), sigma2 = c(1, 1),
                    D = NULL, df = Inf, model = "intslope", theta0 = -3, theta1 = 1,
                    censoring = TRUE, censlam = exp(-3), truncation = TRUE,
                    trunctime = (ntms - 1) + 0.1) {

  K <- nrow(beta)
  if ((K != length(gamma.y)) || (K != length(sigma2))) {
    stop("Incompatible dimensions: number of outcomes\n")
  }

  if (K == 1) {
    stop("Error: this function on simulates multivariate data")
  }

  if (model == "intslope") {
    rk <- 2
  } else if (model == "int") {
    rk <- 1
  } else {
    stop(paste("Unknown model:", model))
  }
  r <- rk * K
  pk <- ncol(beta)
  p <- pk * K
  q <- ifelse(is.null(gamma.x), 0, length(gamma.x))

  if (missing(D) | is.null(D)) {
    D <- diag(r)
  }

  if (length(D) != r^2) {
    stop("Incompatible dimensions: covariance matrix\n")
  }

  if (!isSymmetric(D)) {
    stop("Covariance matrix is not symmetric")
  }

  if (any(eigen(D)$values < 0) || (det(D) <= 0)) {
    stop("Covariance matrix must be positive semi-definite")
  }

  ctsx <- rnorm(n)
  #binx <- runif(n, -sqrt(3), sqrt(3))
  binx <- rbinom(n, 1, 0.5)
  V <- cbind(ctsx, binx)
  id <- 1:n
  idl <- rep(id, each = ntms)
  ctsxl <- rep(ctsx, each = ntms)
  binxl <- rep(binx, each = ntms)
  time <- rep(0:(ntms-1), length = n*ntms)
  X <- cbind(intercept = 1, ltime = time, ctsxl, binxl) # common to all K-variates
  if (df == Inf) {
    b <- MASS::mvrnorm(n, mu = rep(0, r), Sigma = D)
  } else {
    b <- mvtnorm::rmvt(n, sigma = D, df = df)
  }
  bl <- b[rep(1:n, each = ntms), ]
  Z <- matrix(0, nrow =  rk, ncol =  length(time))
  for(i in 1:rk) {
    Z[i, ] <- time^(i-1)
  }
  Y <- matrix(NA, nrow =  n*ntms, ncol =  K)
  for(k in 1:K) {
    Zb <- t(Z) * bl[, ((k-1)*rk + 1):(k*rk), drop = FALSE]
    Y[, k] <- (X %*% beta[k, ]) + rowSums(Zb) + sqrt(sigma2[k]) * rnorm(n*ntms)
  }
  b0k <- b[, ((0:(K-1))*rk + 1)]
  if (model == "intslope") {
    b1k <- b[, ((1:K)*rk)]
  } else {
    b1k <- matrix(0, nrow =  n, ncol =  K)
  }
  Vgam <- V %*% gamma.x
  cens <- rep(1, n)

  U <- runif(n)
  if (model == "int") {
    survtime <- -log(U) / exp(theta0 + Vgam + b0k %*% gamma.y)
  } else {
    ii <- ((theta1 + b1k %*% gamma.y) < 0) &
      (U < exp(exp(theta0 + Vgam + b0k %*% gamma.y) / (theta1 + b1k %*% gamma.y)))
    survtime <- rep(0, n)
    survtime[ii] <- Inf
    survtime[!ii] <- log(
      1 - (theta1 + b1k[!ii, ] %*% gamma.y) * log(U[!ii]) / exp(
        theta0 + Vgam[!ii] + b0k[!ii, ] %*% gamma.y)) / (theta1 + b1k[!ii, ] %*% gamma.y)
  }

  if (censoring) {
    censtime <- -log(runif(n)) / censlam
  } else {
    censtime <- rep(Inf, n)
  }

  if (truncation) {
    censtime <- pmin(censtime, trunctime)
  }

  ii <- (censtime < survtime)
  survtime[ii] <- censtime[ii]
  cens[ii] <- 0
  cens[1] <- 1 # EM algorithm currently only works with this
  ls <- rep(survtime, each = ntms)
  Y <- Y[ls > time, ]
  X <- X[ls > time, ]
  idl <- idl[ls > time]
  time <- time[ls > time]
  cat(paste0(round(100 * sum(cens) / n, 1), "% experienced event\n"))

  list(longdat = data.frame(id = idl, Y = Y, time, X),
       survdat = data.frame(id, survtime, cens, V))

}





