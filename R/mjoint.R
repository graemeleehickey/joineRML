#' Fit a joint model to time-to-event data and multivariate longitudinal data
#'
#' @description This function fits the joint model proposed by Henderson et al.
#'   (2000), but extended to the case of multiple continuous longitudinal
#'   measures. The time-to-event data is modelled using a Cox proportional
#'   hazards regression model with time-varying covariates. The multiple
#'   longitudinal outcomes are modelled using a multivariate version of the
#'   Laird and Ware linear mixed model. The association is captured by a
#'   multivariate latent Gaussian process. The model is estimated using a
#'   (quasi-) Monte Carlo Expectation Maximization (MCEM) algorithm.
#'
#' @param formLongFixed a list of formulae for the fixed effects component of
#'   each longitudinal outcome. The left hand-hand side defines the response,
#'   and the right-hand side specifies the fixed effect terms. If a single
#'   formula is given (either as a list of length 1 or a formula), then it is
#'   assumed that a standard univariate joint model is being fitted.
#' @param formLongRandom a list of one-sided formulae specifying the model for
#'   the random effects effects of each longitudinal outcome. The length of the
#'   list must be equal to \code{formLongFixed}.
#' @param formSurv a formula specifying the proportional hazards regression
#'   model (not including the latent association structure). See
#'   \code{\link[survival]{coxph}} for examples.
#' @param data a list of \code{data.frame} objects for each longitudinal outcome
#'   in which to interpret the variables named in the \code{formLongFixed} and
#'   \code{formLongRandom}. The \code{list} structure enables one to include
#'   multiple longitudinal outcomes with different measurement protocols. If the
#'   multiple longitudinal outcomes are measured at the same time points for
#'   each patient, then a \code{data.frame} object can be given instead of a
#'   \code{list}. It is assumed that each data frame is in long format. \code{tibble}
#'   objects are automatically converted to plain \code{data.frame} objects.
#' @param survData a \code{data.frame} in which to interpret the variables named
#'   in the \code{formSurv}. This is optional, and if not given, the required
#'   data is searched for in \code{data}. Default is \code{survData=NULL}.
#' @param timeVar a character string indicating the time variable in the linear
#'   mixed effects model. If there are multiple longitudinal outcomes and the
#'   time variable is labelled differently in each model, then a character
#'   string vector can be given instead.
#' @param inits a list of initial values for some or all of the parameters
#'   estimated in the model. Default is \code{NULL}, with initial values
#'   estimated using separate multivariate linear mixed effects and Cox
#'   proportional hazard regression models.
#' @param verbose logical: if \code{TRUE}, the parameter estimates and other
#'   convergence statistics are value are printed at each iteration of the MCEM
#'   algorithm. Default is \code{FALSE}.
#' @param pfs logical: if \code{TRUE}, then assuming the MCEM algorithm has
#'   converged, post-fit statistics including the posterior means and variances
#'   of the random effects, and the approximate standard errors are calculated
#'   and returned as part of the model object. Default is \code{TRUE}. If
#'   \code{FALSE}, then these additional calculations are not performed, which
#'   can reduce the overall computational time. This option is intended to be
#'   used with computationally intensive routines such as simulation and
#'   bootstrap standard error estimation where these calculations are not
#'   required.
#' @param control a list of control values with components: \describe{
#'
#'   \item{\code{nMC}}{integer: the initial number of Monte Carlo samples to be
#'   used for integration in the burn-in phase of the MCEM. Default is
#'   \code{nMC=}100\emph{K}.}
#'
#'   \item{\code{nMCscale}}{integer: the scale factor for the increase in Monte
#'   Carlo size when Monte Carlo has not reduced from the previous iteration.
#'   Default is \code{nMCscale=5}.}
#'
#'   \item{\code{nMCmax}}{integer: the maximum number of Monte Carlo samples
#'   that the algorithm is allowed to reach. Default is \code{nMCmax=20000}.}
#'
#'   \item{\code{burnin}}{integer: the number of iterations for 'burn-in' phase
#'   of the optimization algorithm. It is computationally inefficient to use a
#'   large number of Monte Carlo samples early on until one is approximately
#'   near the maximum likelihood estimate. Default is \code{burnin=}100\emph{K}
#'   for \code{type='antithetic'} or \code{type='montecarlo'} and
#'   \code{burnin=}5 for \code{type='sobol'} or \code{type='halton'}. For
#'   standard methods, such a large burn-in will generally be unnecessary and
#'   can be reduced on an application-specific basis.}
#'
#'   \item{\code{mcmaxIter}}{integer: the maximum number of MCEM algorithm
#'   iterations allowed. Default is \code{mcmaxIter=burnin+200}.}
#'
#'   \item{\code{convCrit}}{character string: the convergence criterion to be
#'   used. See \strong{Details}.}
#'
#'   \item{\code{gammaOpt}}{character string: by default (\code{gammaOpt='NR'}),
#'   \eqn{\gamma} is updated using a one-step Newton-Raphson iteration, with the
#'   Hessian matrix calculated exactly. If \code{gammaOpt='GN'}, a Gauss-Newton
#'   algorithm-type iteration is implemented, where the Hessian matrix is
#'   approximated based on calculations similar to those used for calculating
#'   the empirical information matrix? If it is used, then the step-length is
#'   adjusted by a nominal scaling parameter of 0.5 in order to reduce the
#'   chance of over-shooting the maximizer.}
#'
#'   \item{\code{tol0}}{numeric: tolerance value for convergence in the
#'   parameters; see \strong{Details}. Default is \code{tol0=1e-03}.}
#'
#'   \item{\code{tol1}}{numeric: tolerance value for convergence in the
#'   parameters; see \strong{Details}. Default is \code{tol1=1e-03}.}
#'
#'   \item{\code{tol2}}{numeric: tolerance value for convergence in the
#'   parameters; see \strong{Details}. Default is \code{tol2=5e-03} for
#'   \code{type='antithetic'} or \code{type='montecarlo'} and \code{tol2=1e-03}
#'   for \code{type='sobol'} or \code{type='halton'}.}
#'
#'   \item{\code{tol.em}}{numeric: tolerance value for convergence in the
#'   multivariate linear mixed model (MV-LMM). When \eqn{K > 1}, the optimal
#'   initial parameters are those from the MV-LMM, which is estimated using a
#'   separate EM algorithm. Since both the E- and M-steps are available in
#'   closed-form, this algorithm convergences relatively rapidly with a high
#'   precision. Default is min(\code{1e-04}, \code{tol2}).}
#'
#'   \item{\code{rav}}{numeric: threshold when using \code{convCrit='sas'} that
#'   applies absolute change (when \eqn{<}\code{rav}) or relative change (when
#'   \eqn{\ge}\code{rav}) criterion; see \strong{Details}. Default is
#'   \code{0.1}, which is an order of magnitude higher than the SAS
#'   implementation.}
#'
#'   \item{\code{type}}{character: type of Monte Carlo integration method to
#'   use. Options are \describe{
#'
#'   \item{\code{type='montecarlo'}}{Vanilla Monte Carlo sampling.}
#'
#'   \item{\code{type='antithetic'}}{Variance reduction method using antithetic
#'   simulation. This is the default option.}
#'
#'   \item{\code{type='sobol'}}{Quasi-Monte Carlo with a low
#'   deterministic Sobol sequence with Owen-type scrambling.}
#'
#'   \item{\code{type='halton'}}{Quasi-Monte Carlo with a low deterministic
#'   Halton sequence.}
#'
#'   }}
#'
#'   }
#' @param ... options passed to the \code{control} argument.
#'
#' @details Function \code{mjoint} fits joint models for time-to-event and
#'   multivariate longitudinal data. A review of relevant statistical
#'   methodology for joint models of multivariate data is given in Hickey et al.
#'   (2016). This is a direct extension of the models developed in the seminal
#'   works of Wulfsohn and Tsiatis (1997) and Henderson et al. (2000), with the
#'   calculations based largely on Lin et al. (2002) who also extended the model
#'   to multivariate joint data. A more detailed explanation of the model
#'   formulation is given in the technical vignette. Each longitudinal outcome
#'   is modelled according to a linear mixed model (LMM), akin to the models fit
#'   by \code{\link[nlme]{lme}}, with independent and identically distributed
#'   Gaussian errors. The latent term in each model (specified by
#'   \code{formLongRandom}) is a linear combination of subject-specific
#'   zero-mean Gaussian random effects with outcome-specific variance
#'   components. We denote these as \eqn{W_{i1}(t, b_{i1})}, \eqn{W_{i2}(t,
#'   b_{i2})}, \eqn{\dots}, \eqn{W_{iK}(t, b_{iK})}, for the \emph{K}-outcomes.
#'   Usually, \eqn{W(t, b)} will be specified as either \eqn{b_0} (a
#'   random-intercepts model) or \eqn{b_0 + b_1t} (a random-intercepts and
#'   random-slopes model); however, more general structures are allowed The
#'   time-to-event model is modelled as per the usual Cox model formulation,
#'   with an additional (possibly) time-varying term given by
#'
#'   \deqn{\gamma_{y1} W_{i1}(t, b_{i1}) + \gamma_{y2} W_{i2}(t, b_{i2}) + \dots
#'   + \gamma_{yK} W_{iK}(t, b_{iK}),}
#'
#'   where \eqn{\gamma_y} is a parameter vector of proportional latent
#'   association parameters of length \emph{K} for estimation.
#'
#'   The optimization routine is based on a Monte Carlo Expectation Maximization
#'   algorithm (MCEM) algorithm, as described by Wei and Tanner (1990). With
#'   options for using antithetic simulation for variance reduction in the Monte
#'   Carlo integration, or quasi-Monte Carlo based on low order deterministic
#'   sequences.
#'
#' @section Convergence criteria:
#'
#'   The routine internally scales and centers data to avoid overflow in the
#'   argument to the exponential function. These actions do not change the
#'   result, but lead to more numerical stability. Several convergence criteria
#'   are available: \describe{
#'
#'   \item{\code{abs}}{the maximum absolute parameter change is
#'   \eqn{<}\code{tol0}. The baseline hazard parameters are not included in this
#'   convergence statistic.}
#'
#'   \item{\code{rel}}{the maximum (absolute) relative parameter change is
#'   \eqn{<}\code{tol2}. A small value (\code{tol1}) is added to the denominator
#'   of the relative change statistic to avoid numerical problems when the
#'   parameters are close to zero.}
#'
#'   \item{\code{either}}{\emph{either} the \code{abs} or \code{rel} criteria
#'   are satisfied.}
#'
#'   \item{\code{sas}}{if \eqn{| \theta_p | < }\code{rav}, then the \code{abs}
#'   criteria is applied for the \emph{l}-th parameter; otherwise, \code{rel} is
#'   applied. This is the approach used in the SAS EM algorithm program for
#'   missing data.}
#'
#'   }
#'
#'   Due to the Monte Caro error, the algorithm could spuriously declare
#'   convergence. Therefore, we require convergence to be satisfied for 3
#'   consecutive iterations. The algorithm starts with a low number of Monte
#'   Carlo samples in the 'burn-in' phase, as it would be computationally
#'   inefficient to use a large sample whilst far away from the true maximizer.
#'   After the algorithm moves out of this adaptive phase, it uses an automated
#'   criterion based on the coefficient of variation of the relative parameter
#'   change of the last 3 iterations to decide whether to increase the Monte
#'   Carlo sample size. See the technical vignette and Ripatti et al. (2002) for
#'   further details.
#'
#' @section Standard error estimation:
#'
#'   Approximate standard errors (SEs) are calculated (if \code{pfs=TRUE}).
#'   These are based on the empirical observed information function (McLachlan &
#'   Krishnan, 2008). Through simulation studies, we have found that this
#'   approximation does not work particularly well for \eqn{n < 100} (where
#'   \emph{n} is the number of subjects). In these cases, one would need to
#'   appeal to the bootstrap SE estimation approach. However, in practice, the
#'   reliability of the approximate SEs will depend of a multitude of factors,
#'   including but not limited to, the average number of repeated measurements
#'   per subject, the total number of events, and the convergence of the MCEM
#'   algorithm.
#'
#'   Bootstrap SEs are also available, however they are not calculated using the
#'   \code{mjoint} function due to the intense computational time. Instead, a
#'   separate function is available: \code{bootSE}, which takes the fitted joint
#'   model as its main argument. Given a fitted joint model (of class
#'   \code{mjoint}) and a bootstrap fit object (of class \code{bootSE}), the SEs
#'   reported in the model can be updated by running
#'   \code{summary(fit_obj,boot_obj)}. For details, consult the
#'   \code{\link{bootSE}} documentation.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords multivariate survival methods
#' @seealso \code{\link{mjoint.object}}, \code{\link{bootSE}},
#'   \code{\link{plot.mjoint}}, \code{\link{summary.mjoint}},
#'   \code{\link{getVarCov.mjoint}}, \code{\link{simData}}.
#'
#' @references
#'
#' Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
#' measurements and event time data. \emph{Biostatistics.} 2000; \strong{1(4)}:
#' 465-480.
#'
#' Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. Joint modelling of
#' time-to-event and multivariate longitudinal outcomes: recent developments and
#' issues. \emph{BMC Med Res Methodol.} 2016; \strong{16(1)}: 117.
#'
#' Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the joint
#' analysis of time-to-event and multiple longitudinal variables. \emph{Stat
#' Med.} 2002; \strong{21}: 2369-2382.
#'
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second
#' Edition. Wiley-Interscience; 2008.
#'
#' Ripatti S, Larsen K, Palmgren J. Maximum likelihood inference for
#' multivariate frailty models using an automated Monte Carlo EM algorithm.
#' \emph{Lifetime Data Anal.} 2002; \strong{8}: 349-360.
#'
#' Wei GC, Tanner MA. A Monte Carlo implementation of the EM algorithm and the
#' poor man's data augmentation algorithms. \emph{J Am Stat Assoc.} 1990;
#' \strong{85(411)}: 699-704.
#'
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
#' measured with error. \emph{Biometrics.} 1997; \strong{53(1)}: 330-339.
#'
#' @import stats survival
#'
#' @return An object of class \code{mjoint}. See \code{\link{mjoint.object}} for
#'   details.
#' @export
#'
#' @examples
#' # Fit a classical univariate joint model with a single longitudinal outcome
#' # and a single time-to-event outcome
#'
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#'
#' set.seed(1)
#' fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
#'     formLongRandom = ~ time | num,
#'     formSurv = Surv(fuyrs, status) ~ age,
#'     data = hvd,
#'     timeVar = "time",
#'     control = list(nMCscale = 2, burnin = 5)) # controls for illustration only
#' summary(fit1)
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
#' summary(fit2)
#' }
#'
#' \dontrun{
#' # Fit a univariate joint model and compare to the joineR package
#'
#' data(pbc2)
#' pbc2$log.b <- log(pbc2$serBilir)
#'
#' # joineRML package
#' fit.joineRML <- mjoint(
#'     formLongFixed = list("log.bil" = log.b ~ year),
#'     formLongRandom = list("log.bil" = ~ 1 | id),
#'     formSurv = Surv(years, status2) ~ age,
#'     data = pbc2,
#'     timeVar = "year")
#' summary(fit.joineRML)
#'
#' # joineR package
#' pbc.surv <- UniqueVariables(pbc2, var.col = c("years", "status2"),
#'                             id.col = "id")
#' pbc.long <- pbc2[, c("id", "year", "log.b")]
#' pbc.cov <- UniqueVariables(pbc2, "age", id.col = "id")
#' pbc.jd <- jointdata(longitudinal = pbc.long, baseline = pbc.cov,
#'                     survival = pbc.surv, id.col = "id", time.col = "year")
#' fit.joineR <- joint(data = pbc.jd,
#'     long.formula = log.b ~ 1 + year,
#'     surv.formula = Surv(years, status2) ~ age,
#'     model = "int")
#' summary(fit.joineR)
#' }
mjoint <- function(formLongFixed, formLongRandom, formSurv, data, survData = NULL,
                   timeVar, inits = NULL, verbose = FALSE, pfs = TRUE,
                   control = list(), ...) {

  #*****************************************************
  # Preamble
  #*****************************************************

  time.start <- Sys.time()
  Call <- match.call()
  balanced <- FALSE # assume unless proven o/w

  # Formulas do not need to be given as lists if K=1
  if (!is.list(formLongFixed)) {
    balanced <- TRUE
    formLongFixed <- list(formLongFixed)
    formLongRandom <- list(formLongRandom)
    K <- 1
  } else {
    K <- length(formLongFixed)
  }

  # Data does not need to a list if K=1
  # if K>1 and not a list, assume data balanced
  if (!("list" %in% class(data))) {
    balanced <- TRUE
    data <- list(data)
    if (K > 1) {
      for (k in 2:K) {
        data[[k]] <- data[[1]]
      }
    }
  } else {
    balanced <- (length(unique(data)) == 1)
  }
  if (length(data) != K) {
    stop(paste("The number of datasets expected is K =", K))
  }

  # Strip `tibble` class from data (breaks indexing by column) if present
  data <- lapply(data, function(d) {
    if (any(c("tbl_df", "tbl") %in% class(d))) {
      return(as.data.frame(d))
    } else {
      return(d)
    }
  })

  id <- as.character(nlme::splitFormula(formLongRandom[[1]], "|")[[2]])[2]
  n <- length(unique(data[[1]][, id]))

  # Incase timeVar not a vector when K>1
  if (length(timeVar) == 1 & (K > 1)) {
    timeVar <- rep(timeVar, K)
  }
  if (length(timeVar) != K) {
    stop(paste("The length of timeVar must equal", K))
  }

  # Order the data + id -> factor
  for (k in 1:K) {
    data[[k]] <- data[[k]][order(xtfrm(data[[k]][, id]), data[[k]][, timeVar[k]]), ]
    data[[k]][, id] <- as.factor(data[[k]][, id])
    data[[k]] <- droplevels(data[[k]])
  }

  # Check the same patients measured at least once for each marker
  if (K > 1) {
    uniq.ids <- list(sort(unique(data[[1]][, id])))
    for (k in 2:K) {
      uniq.ids[[k]] <- sort(unique(data[[k]][, id]))
      if (!identical(uniq.ids[[k - 1]], uniq.ids[[k]])) {
        stop("Every subject must have at least one measurement per each outcome")
      }
    }
  }

  #*****************************************************
  # Control parameters
  #*****************************************************

  con <- list(nMC = 100*K, nMCscale = 5, nMCmax = 20000, burnin = 100*K,
              mcmaxIter = 100*K + 200, convCrit = "sas", gammaOpt = "NR",
              tol0 = 1e-03, tol1 = 1e-03, tol2 = 5e-03, tol.em = 1e-04,
              rav = 0.1, type = "antithetic")
  nc <- names(con)
  control <- c(control, list(...))
  con[(conArgs <- names(control))] <- control

  # Checks and over-rides

  # Default for tol.em
  if (("tol.em" %in% names(control)) || ("tol2" %in% names(control))) {
    con$tol.em <- min(con$tol.em, con$tol2)
  }
  # Don't allow algorithm to terminate before end of burnin
  if (("burnin" %in% names(control)) && !("mcmaxIter" %in% names(control))) {
    con$mcmaxIter <- con$burnin + 200
  }
  # Separate defaults for QMC methods
  if (con$type %in% c("sobol", "halton")) {
    # reduce burn-in
    if (!("burnin" %in% names(control))) {
      con$burnin <- 5
    }
    # reduce tol2
    if (!("tol2" %in% names(control))) {
      con$tol2 <- 1e-03
    }
  }
  # If any unmatched arguments given
  if (length(unmatched <- conArgs[!(conArgs %in% nc)]) > 0) {
    warning("Unknown arguments passed to 'control': ", paste(unmatched, collapse = ", "))
  }

  #*****************************************************
  # Multivariate longitudinal data
  #*****************************************************

  lfit <- list()
  mf.fixed <- list()
  yik <- list()
  Xik <- list()
  nk <- vector(length = K)
  Xik.list <- list()
  nik.list <- list()
  Zik <- list()
  Zik.list <- list()

  for (k in 1:K) {

    # List of K seperate longitudinal model fits
    lfit[[k]] <- nlme::lme(fixed = formLongFixed[[k]], random = formLongRandom[[k]],
                           data = data[[k]], method = "ML",
                           control = nlme::lmeControl(opt = "optim"))
    lfit[[k]]$call$fixed <- eval(lfit[[k]]$call$fixed)

    # Model frames
    mf.fixed[[k]] <- model.frame(lfit[[k]]$terms,
                                 data[[k]][, all.vars(formLongFixed[[k]])])

    # Longitudinal outcomes
    yik[[k]] <- by(model.response(mf.fixed[[k]], "numeric"), data[[k]][, id], as.vector)

    # X design matrix
    Xik[[k]] <- data.frame("id" = data[[k]][, id],
                           model.matrix(formLongFixed[[k]], data[[k]]))

    # n_k (number of obs per each k)
    nk[k] <- nrow(Xik[[k]])

    # X design matrix (list)
    Xik.list[[k]] <- by(Xik[[k]], Xik[[k]]$id, function(u) {
      as.matrix(u[, -1])
    })

    # Sample size for each subject
    nik.list[[k]] <- by(Xik[[k]], Xik[[k]]$id, nrow)

    # Z design matrix
    ffk <- nlme::splitFormula(formLongRandom[[k]], "|")[[1]]
    Zik[[k]] <- data.frame("id" = data[[k]][, id], model.matrix(ffk, data[[k]]))

    # Z design matrix (list)
    Zik.list[[k]] <- by(Zik[[k]], Zik[[k]]$id, function(u) {
      as.matrix(u[, -1])
    })

  }

  # Flatten lists to length = n
  yi <- sapply(names(yik[[1]]), function(i) {
    unlist(lapply(yik, "[[", i))
  },
  USE.NAMES = TRUE, simplify = FALSE)

  Xi <- sapply(names(Xik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Xik.list, "[[", i)))
  },
  USE.NAMES = TRUE, simplify = FALSE)

  Zi <- sapply(names(Zik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Zik.list, "[[", i)))
  },
  USE.NAMES = TRUE, simplify = FALSE)

  Zit <- lapply(Zi, t)

  # t(X) %*% X [summed over i and inverted]
  XtXi <- lapply(Xi, crossprod)
  XtX.inv <- solve(Reduce("+", XtXi))

  # t(X) %*% y [for each i]
  Xtyi <- mapply(function(x, y) {
    crossprod(x, y)
  },
  x = Xi, y = yi,
  SIMPLIFY = FALSE)

  # t(X) %*% Z [for each i]
  XtZi <- mapply(function(x, z) {
    crossprod(x, z)
  },
  x = Xi, z = Zi,
  SIMPLIFY = FALSE)

  nik <- sapply(names(nik.list[[1]]), function(i) {
    unlist(lapply(nik.list, "[[", i))
  },
  USE.NAMES = TRUE, simplify = FALSE)

  # Number of fixed and random effects
  p <- sapply(1:K, function(i) {
    ncol(Xik[[i]]) - 1
  })

  r <- sapply(1:K, function(i) {
    ncol(Zik[[i]]) - 1
  })

  l <- list(yi = yi, Xi = Xi, Zi = Zi, Zit = Zit, nik = nik, yik = yik,
            Xik.list = Xik.list, Zik.list = Zik.list, XtX.inv = XtX.inv,
            Xtyi = Xtyi, XtZi = XtZi, p = p, r = r, K = K, n = n, nk = nk)

  #*****************************************************
  # Time-to-event data
  #*****************************************************

  if (is.null(survData)) {
    survData <- data[[1]][!duplicated(data[[1]][, id]), ]
  }
  survdat <- survData

  sfit <- survival::coxph(formSurv, data = survdat, x = TRUE)
  q <- ncol(sfit$x)
  sfit.start <- survival::survfit(sfit)
  tj <- sfit.start$time[sfit.start$n.event > 0]
  nev <- sfit.start$n.event[sfit.start$n.event > 0]
  nev.uniq <- length(tj)

  survdat2 <- data.frame(survdat[, id], sfit$x, sfit$y[, 1], sfit$y[, 2])
  xcenter <- NULL
  if (q > 0) {
    xcenter <- apply(survdat2[2:(q + 1)], 2, mean)
    survdat2[2:(q + 1)] <- scale(survdat2[2:(q + 1)],
                                 center = xcenter, scale = FALSE)
  }
  colnames(survdat2)[c(1, (q + 2):(q + 3))] <- c("id", "T", "delta")
  survdat2$tj.ind <- sapply(1:n, function(i) {
    sum(tj <= survdat2$T[i])
  })
  survdat2.list <- by(survdat2, survdat2$id, list)
  if (q > 0) {
    V <- by(survdat2, survdat2$id, function(u) {
      unlist(u[, 2:(q + 1)])
    })
  } else {
    V <- by(survdat2, survdat2$id, function(u) {
      0
    })
  }

  # Collect together as inputs for EM algorithm
  t <- list(V = V, survdat2 = survdat2, survdat2.list = survdat2.list,
            q = q, nev = nev, nev.uniq = nev.uniq, xcenter = xcenter,
            tj = tj)

  # Longitudinal data should not be recorded *after* event time
  for (k in 1:K) {
    for (i in survdat2$id) {
      if (max(data[[k]][data[[k]][, id] == i, timeVar[k]]) > survdat2[survdat2$id == i, "T"]) {
        stop("Longitudinal measurements should not be recorded after the event time")
      }
    }
  }

  # Separate models log-likelihood
  # NB: adjusts for number of events as sfit loglik is a partial estimate
  log.lik0 <- sum(sapply(lfit, logLik)) + sfit$loglik[ifelse(q > 0, 2, 1)] - sfit$nevent
  log.lik <- log.lik0

  #*****************************************************
  # Z(t_j) for unique failure times t_j
  #*****************************************************

  # Subjects who have a censoring time before the first event time
  # must be handled carefully... they don't fully contribute to the
  # likelihood. We allow them to have 1 row in the dataset to ensure
  # non-empty lists, but effectively discount this data in the EM-algorithm
  Zdat.fail <- data.frame(
    "id" = rep(unique(survdat2$id), pmax(survdat2$tj.ind, 1)),
    "time" = unlist(sapply(1:n, function(i) {
      tj[1:max(survdat2$tj.ind[i], 1)]
    },
    simplify = FALSE))
  )
  names(Zdat.fail)[1] <- id

  Zik.fail <- list()
  Zik.fail.list <- list()

  for(k in 1:K) {

    if (ncol(Zdat.fail) == 2) {
      names(Zdat.fail)[2] <- timeVar[k]
    }

    # Z design matrix
    ffk <- nlme::splitFormula(formLongRandom[[k]], "|")[[1]]
    Zik.fail[[k]] <- data.frame("id" = Zdat.fail[, id], model.matrix(ffk, Zdat.fail))

    # Z design matrix (list)
    Zik.fail.list[[k]] <- by(Zik.fail[[k]], Zik.fail[[k]]$id, function(u) {
      as.matrix(u[, -1])
    })

  }

  # Flatten list to length = n
  Zi.fail <- sapply(names(Zik.fail.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Zik.fail.list, "[[", i)))
  },
  USE.NAMES = TRUE, simplify = FALSE)
  Zit.fail <- lapply(Zi.fail, t) # transpose

  # Side-stacked identity matrices for calculating Zb
  IW.fail <- by(survdat2, survdat2$id, function(u) {
    do.call("cbind", lapply(1:K, function(i) diag(max(u$tj.ind, 1))))
  })

  # Collect together as inputs for EM algorithm
  z <- list(Zi.fail = Zi.fail, Zit.fail = Zit.fail, IW.fail = IW.fail)

  #*****************************************************
  # Initial values
  #*****************************************************

  if (!is.null(inits)) {
    if (!is.list(inits)) {
      stop("inits must be a list.\n")
    }
    theta.names <- c("D", "beta", "sigma2", "gamma", "haz")
    if (length(unmatched <- names(inits)[!(names(inits) %in% theta.names)]) > 0) {
      warning("Unknown initial parameters passed to 'inits': ", paste(unmatched, collapse = ", "))
    }
  }

  inits.long <- initsLong(lfit = lfit, inits = inits, l = l, z = z, K = K, p = p, r = r,
                          tol.em = con$tol.em, verbose = verbose)

  inits.surv <- initsSurv(data = data, lfit = lfit, sfit = sfit, survdat2 = survdat2,
                          formSurv = formSurv, id = id, timeVar = timeVar,
                          K = K, q = q, balanced = balanced, inits = inits)

  D <- inits.long$D
  beta <- inits.long$beta
  sigma2 <- inits.long$sigma2
  gamma <- inits.surv$gamma
  haz <- inits.surv$haz

  theta <- list("D" = D, "beta" = beta, "sigma2" = sigma2,
                "gamma" = gamma, "haz" = haz)

  #*****************************************************
  # Run EM algorithm
  #*****************************************************

  rm(yi, Xi, Zi, Zit, yik, Xik.list, Zik.list, XtX.inv, Xtyi, XtZi,
     V, survdat2, survdat2.list, nev, nev.uniq, xcenter, tj,
     Zdat.fail, Zi.fail, Zit.fail, IW.fail)

  all.iters <- list()
  conv.track <- rep(FALSE, con$mcmaxIter)
  Delta.vec <- rep(NA, con$mcmaxIter)
  ll.hx <- rep(NA, con$mcmaxIter)
  cv.old <- 0
  time.em.start <- Sys.time()

  nMC <- con$nMC
  nmc.iters <- c()

  for (it in 1:(con$mcmaxIter)) {

    if (verbose) {
      cat("-------------------------------------------------------------\n\n")
      cat(paste0("Iteration: ", it, "\n\n"))
    }

    nmc.iters <- c(nmc.iters, nMC)

    stepUpdate <- stepEM(theta = theta, l = l, t = t, z = z,
                         nMC = nMC, verbose = verbose, gammaOpt = con$gammaOpt,
                         pfs = FALSE, type = con$type)
    theta.new <- stepUpdate$theta.new
    log.lik.new <- stepUpdate$ll
    ll.hx[it] <- log.lik.new

    all.iters[[it]] <- theta.new
    if (verbose) {
      print(theta.new[-which(names(theta.new) == "haz")])
      #print(theta.new)
    }

    # Stepwise convergence
    conv.status <- convMonitor(theta = theta, theta.new = theta.new,
                               log.lik = log.lik, log.lik.new = log.lik.new,
                               con = con, verbose = verbose)
    conv.track[it] <- conv.status$conv
    Delta.vec[it] <- conv.status$max.reldelta.pars
    log.lik <- log.lik.new

    if (it >= con$burnin + 3) {
      # require convergence condition to be satisfied 3 iterations
      # in a row + cannot converge during burn-in stage
      conv <- all(conv.track[(it - 2):it])
    } else {
      conv <- FALSE
    }

    if (!conv && (it >= con$burnin)) {
      # Ripatti decision-rule for nMC increase using CV statistics
      cv <- ifelse(it >= 3, sd(Delta.vec[(it - 2):it]) / mean(Delta.vec[(it - 2):it]), 0)
      if (verbose) {
        cat(paste("CV statistic (old) =", round(cv.old, 6), "\n"))
        cat(paste("CV statistic (new) =", round(cv, 6), "\n\n"))
      }
      ripatti <- (cv > cv.old)
      cv.old <- cv
      # Check for drift in log-likelihood with using QMC
      lldrift <- (con$type %in% c("sobol", "halton")) && (conv.status$rel.ll < 0)
      if (ripatti || lldrift) {
        nMC.old <- nMC
        nMC <- min(nMC + floor(nMC / con$nMCscale), con$nMCmax)
      }
      if (verbose) {
        cat(paste("Using number of MC nodes:", nMC, "\n\n"))
      }
    }

    # Once converged: calculate SEs and posterior REs (means + variances)
    if (conv) {
      message("EM algorithm has converged!\n")
      theta <- theta.new
      if (pfs) {
        message("Calculating post model fit statistics...\n")
        postFitCalcs <- stepEM(theta = theta, l = l, t = t, z = z,
                               nMC = nMC, verbose = FALSE, gammaOpt = "NR",
                               pfs = TRUE, type = con$type)
      }
      break
    } else {
      theta <- theta.new
    }
    if ((it == con$mcmaxIter) & !conv) {
      message("Failed to converge!")
    }
    utils::flush.console()

  }

  # History
  hx.beta <- sapply(all.iters, function(x) x$beta)
  rownames(hx.beta) <- paste0("long_", rownames(hx.beta))
  hx.gamma <- sapply(all.iters, function(x) x$gamma)
  if ((q + K) == 1) {
    hx.gamma <- matrix(hx.gamma, nrow = 1)
  }
  rownames(hx.gamma) <- paste0("surv_", rownames(hx.gamma))
  hx.D <- sapply(all.iters, function(x) x$D[lower.tri(x$D, diag = TRUE)])
  if (sum(r) == 1) {
    hx.D <- matrix(hx.D, nrow = 1)
  }
  ltri <- lower.tri(theta$D, diag = TRUE)
  rownames(hx.D) <- paste0("D_", row = row(theta$D)[ltri], ",",
                           col = col(theta$D)[ltri])
  hx.haz <- sapply(all.iters, function(x) x$haz)
  rownames(hx.haz) <- paste0("haz_", 1:nrow(hx.haz))
  hx.sigma2 <- sapply(all.iters, function(x) x$sigma2)
  if (K == 1) {
    hx.sigma2 <- matrix(hx.sigma2, nrow = 1)
  }
  if (K > 1) {
    rownames(hx.sigma2) <- paste0("sigma2_", 1:K)
  } else {
    rownames(hx.sigma2) <- "sigma2"
  }

  #*****************************************************
  # Output
  #*****************************************************

  out <- list("coefficients" = theta.new)
  out$history <- rbind(hx.beta, hx.gamma, hx.sigma2, hx.D, hx.haz)
  out$nMC.hx <- nmc.iters
  out$formLongFixed <- formLongFixed
  out$formLongRandom <- formLongRandom
  out$formSurv <- formSurv
  out$data <- data          # returns the data in subject-time order
  out$survData <- survData  # returns the data in subject order
  out$timeVar <- timeVar
  out$id <- id
  out$dims <- list("p" = p, "r" = r, "K" = K, "q" = q,
                   "n" = n, "nk" = nk)
  out$sfit <- sfit
  out$lfit <- lfit
  out$log.lik0 <- log.lik0
  out$log.lik <- log.lik
  out$ll.hx <- ll.hx
  out$control <- con
  out$finalnMC <- nMC # not same as control nMC (used for burn-in phase)
  if (conv && pfs) {
    out$Hessian <- postFitCalcs$H
    out$Eb <- postFitCalcs$Eb # Posterior RE means
    out$Vb <- postFitCalcs$Vb # Posterior RE variances
    out$dmats <- list(l = l, t = t, z = z)
  }
  out$call <- Call
  out$conv <- conv

  # Time
  time.end <- Sys.time()
  time.diff <- c("Total" = time.end - time.start,
                 "EM" = time.end - time.em.start)
  out$comp.time <- time.diff
  if (verbose) {
    cat(paste("Total time taken:", round(as.numeric(time.diff[1]), 1),
              attr(time.diff[1], "units"), "\n"))
    cat(paste("EM algorithm took", round(as.numeric(time.diff[2]), 1),
              attr(time.diff[2], "units"), "\n\n"))
  }

  class(out) <- "mjoint"
  invisible(out)

}
