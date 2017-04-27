#' joineRML
#'
#' @description joineRML is an extension of the joineR package for fitting joint
#'   models of time-to-event data and multivariate longitudinal data. The model
#'   fitted in joineRML is an extension of the Wulfsohn and Tsiatis (1997) and
#'   Henderson et al. (2000) models, which is comprised on
#'   \eqn{(K+1)}-sub-models: a Cox proportional hazards regression model (Cox,
#'   1972) and a \emph{K}-variate linear mixed-effects model - a direct
#'   extension of the Laird and Ware (1982) regression model. The model is
#'   fitted using a Monte Carlo Expectation-Maximization (MCEM) algorithm, which
#'   closely follows the methodology presented by Lin et al. (2002).
#'
#' @useDynLib joineRML, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' Cox DR. Regression models and life-tables. \emph{J R Stat Soc Ser B Stat
#' Methodol.} 1972; \strong{34(2)}: 187-220.
#'
#' Laird NM, Ware JH. Random-effects models for longitudinal data.
#' \emph{Biometrics.} 1982; \strong{38(4)}: 963-974.
#'
#' @docType package
#' @name joineRML
NULL
