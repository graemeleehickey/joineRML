#' Fitted \code{mjoint} object
#'
#' @description An object returned by the \code{mjoint} function, inheriting
#'   from class \code{mjoint} and representing a fitted joint model for
#'   multivariate longitudinal and time-to-event data. Objects of this class
#'   have methods for the generic functions \code{coef}, \code{logLik},
#'   \code{plot}, \code{print}, \code{ranef}, \code{fixef}, \code{summary},
#'   \code{AIC}, \code{getVarCov}, \code{vcov}, \code{confint}, \code{sigma},
#'   \code{fitted}, \code{residuals}, and \code{formula}.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords multivariate survival
#' @seealso \code{\link{mjoint}}.
#' @return A list with the following components. \describe{
#'
#'   \item{\code{coefficients}}{a list with the estimated coefficients. The
#'   components of this list are: \describe{
#'
#'   \item{\code{beta}}{the vector of fixed effects for the linear mixed effects
#'   sub-model.}
#'
#'   \item{\code{D}}{the variance-covariance matrix of the random effects.}
#'
#'   \item{\code{sigma2}}{the measurement error standard deviations for the
#'   linear mixed effects sub-model.}
#'
#'   \item{\code{haz}}{the estimated baseline hazard values for each unique
#'   failure time. Note that this is the \emph{centered} hazard, equivalent to
#'   that returned by \code{\link[survival]{coxph.detail}}.}
#'
#'   \item{\code{gamma}}{the vector of baseline covariates for the survival
#'   model and the latent association coefficient parameter estimates.}
#'
#'   }}
#'
#'   \item{\code{history}}{a matrix with parameter estimates at each iteration
#'   of the MCEM algorithm.}
#'
#'   \item{\code{nMC.hx}}{a vector with the number of Monte Carlo samples for
#'   each MCEM algorithm iteration.}
#'
#'   \item{\code{formLongFixed}}{a list of formulae for the fixed effects
#'   component of each longitudinal outcome.}
#'
#'   \item{\code{formLongRandom}}{a list of formulae for the fixed effects
#'   component of each longitudinal outcome. The length of the list will be
#'   equal to \code{formLongFixed}.}
#'
#'   \item{\code{formSurv}}{a formula specifying the proportional hazards
#'   regression model (not including the latent association structure).}
#'
#'   \item{\code{data}}{a list of data.frames for each longitudinal outcome.}
#'
#'   \item{\code{survData}}{a data.frame of the time-to-event dataset.}
#'
#'   \item{\code{timeVar}}{a character string vector of length K denoting the
#'   column name(s) for time in \code{data}.}
#'
#'   \item{\code{id}}{a character string denoting the column name for subject
#'   IDs in \code{data} and \code{survData}.}
#'
#'   \item{\code{dims}}{a list giving the dimensions of model parameters with
#'   components: \describe{
#'
#'   \item{\code{p}}{a vector of the number of fixed effects for each
#'   longitudinal outcome.}
#'
#'   \item{\code{r}}{a vector of the number of random effects for each
#'   longitudinal outcome.}
#'
#'   \item{\code{K}}{an integer of the number of different longitudinal outcome
#'   types.}
#'
#'   \item{\code{q}}{an integer of the number of baseline covariates in the
#'   time-to-event sub-model.}
#'
#'   \item{\code{n}}{an integer of the total number of subjects in the study.}
#'
#'   \item{\code{nk}}{a vector of the number of measurements for each
#'   longitudinal outcome.}
#'
#'   }}
#'
#'   \item{\code{sfit}}{an object of class \code{coxph} for the separate
#'   time-to-event model fit. See \code{\link[survival]{coxph}} for details.}
#'
#'   \item{\code{lfit}}{a list of objects each of class \code{lme} from fitting
#'   separate linear mixed effects models; one per each longitudinal outcome
#'   type. See \code{\link[nlme]{lme}} for details.}
#'
#'   \item{\code{log.lik0}}{the combined log-likelihood from separate sub-model
#'   fits.}
#'
#'   \item{\code{log.lik}}{the log-likelihood from the joint model fit.}
#'
#'   \item{\code{ll.hx}}{a vector of the log-likelihood values for each MCEM
#'   algorithm interaction.}
#'
#'   \item{\code{control}}{a list of control parameters used in the estimation
#'   of the joint model. See \code{\link{mjoint}} for details.}
#'
#'   \item{\code{finalnMC}}{the final number of Monte Carlo samples required
#'   prior to convergence.}
#'
#'   \item{\code{call}}{the matched call.}
#'
#'   \item{\code{conv}}{logical: did the MCEM algorithm converge within the
#'   specified maximum number of iterations?}
#'
#'   \item{\code{comp.time}}{a vector of length 2 with each element an object of
#'   class \code{difftime} that reports the \emph{total} time taken for model
#'   fitting (including all stages) and the time spent in the \emph{EM
#'   algorithm}.}
#'
#'   }
#'
#' @section Post model fit statistics:
#'
#'   If \code{pfs=TRUE}, indicating that post-fit statistics are to be returned,
#'   then the output also includes the following objects. \describe{
#'
#'   \item{\code{vcov}}{the variance-covariance matrix of model parameters, as
#'   approximated by the empirical information matrix, is reported. See
#'   \code{\link{mjoint}} for details.}
#'
#'   \item{\code{SE.approx}}{the square-root of the diagonal of \code{vcov} is
#'   returned, which are estimates of the standard errors for the parameters.}
#'
#'   \item{\code{Eb}}{a matrix with the estimated random effects values for each
#'   subject.}
#'
#'   \item{\code{Vb}}{an array with the estimated variance-covariance matrices
#'   for the random effects values for each subject.}
#'
#'   \item{\code{dmats}}{a list of length 3 containing the design matrices, data
#'   frames, and vectors used in the MCEM algorithm. These are required for
#'   prediction and to calculate the residuals and . The 3 items in the list are
#'   \code{l} (longitudinal data), \code{t} (time-to-event data), and \code{z}
#'   (design matrices expanded over unique failure times). These are not
#'   intended to be extracted by the user.}}
"mjoint.object" <- NULL
