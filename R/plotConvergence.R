#' Plot convergence time series for parameter vectors from an \code{mjoint}
#' object
#'
#' Plot convergence time series for parameter vectors from an \code{mjoint}
#' object.
#'
#' @inheritParams confint.mjoint
#' @param params a string indicating what parameters are to be shown. Options
#'   are \code{params='gamma'} for the time-to-event sub-model covariate
#'   coefficients, including the latent association parameters;
#'   \code{params='beta'} for the longitudinal sub-model fixed effects
#'   coefficients; \code{params='sigma2'} for the residual error variances from
#'   the longitudinal sub-model; \code{params='D'} for the lower triangular
#'   matrix of the variance-covariance matrix of random effects.
#'
#' @references
#'
#' Wei GC, Tanner MA. A Monte Carlo implementation of the EM algorithm and the
#' poor man's data augmentation algorithms. \emph{J Am Stat Assoc.} 1990;
#' \strong{85(411)}: 699-704.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods dplot
#' @seealso \code{\link{plot.mjoint}}, \code{\link[graphics]{plot.default}},
#'   \code{\link[graphics]{par}}, \code{\link[graphics]{abline}}.
#'
#' @importFrom graphics par plot
#' @export
plotConvergence <- function(object, params = "gamma") {

  if (class(object) != "mjoint") {
    stop("Convergence plots require an 'mjoint' model\n")
  }

  his <- object$history
  n.iters <- ncol(his)

  dims <- object$dims
  p <- sum(dims[["p"]])
  q <- dims[["q"]]
  K <- dims[["K"]]
  r <- sum(dims[["r"]])

  # NB: for rows indices of 'his', note number of parameters =
  #     p (beta)
  #     q + K (gamma)
  #     K (sigma2)
  #     r * (r + 1) / 2 (upper D matrix)
  inds <- c(0, cumsum(c(p, q + K, K, r * (r + 1) / 2)))

  if (n.iters == 1) stop("No convergence history found.\n")

  old.par <- par(no.readonly = TRUE)
  nc <- 1

  #--------------------------------------------------------

  # betas
  if (params == "beta") {
    if (p > 3) {nc = ceiling(p / 3)}
    par(mfrow = c(nrow = min(p, 3), ncol = nc))
    for (i in 1:p) {
      plot(his[(inds[1] + 1):(inds[2]), ][i, ], type = "l",
           xlab = "Iteration", ylab = rownames(his[(inds[1] + 1):(inds[2]), ])[i])
    }
  }

  #--------------------------------------------------------

  # gammas
  if (params == "gamma") {
    n.par <- q + K # must always be >= 1
    if (n.par > 3) {nc = ceiling(n.par / 3)}
    par(mfrow = c(nrow = min(n.par, 3), ncol = nc))
    for (i in 1:n.par) {
    plot(his[(inds[2] + 1):(inds[3]), ][i, ], type = "l",
         xlab = "Iteration", ylab = rownames(his[(inds[2] + 1):(inds[3]), ])[i])
    }
  }

  #--------------------------------------------------------

  # sigma2
  if (params == "sigma2") {
    if (K > 3) {nc = ceiling(K / 3)}
    par(mfrow = c(nrow = min(K, 3), ncol = nc))
    for (i in 1:K) {
      plot(his[(inds[3] + 1):(inds[4]), ][i, ], type = "l",
           xlab = "Iteration", ylab = rownames(his[(inds[3] + 1):(inds[4]), ])[i])
    }
  }

  #--------------------------------------------------------

  # D
  if (params == "D") {
    n.par <- r * (r+1) / 2 # upper triangle only
    if (n.par > 3) {nc = ceiling(n.par / 3)}
    par(mfrow = c(nrow = min(n.par, 3), ncol = nc))
    for (i in 1:n.par) {
      plot(his[(inds[4] + 1):(inds[5]), ][i, ], type = "l",
           xlab = "Iteration", ylab = rownames(his[(inds[4] + 1):(inds[5]), ])[i])
    }
  }

  on.exit(par(old.par))

}
