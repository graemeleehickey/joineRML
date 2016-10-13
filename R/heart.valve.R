#' Heart valve surgery
#'
#' This is longitudinal data on an observational study on detecting effects of different heart valves, differing on type of tissue. The data consists of longitudinal measurements on three different heart function outcomes, after surgery occurred. There are several baseline covariates available, and also survival data.
#'
#' @usage data(heart.valve)
#' @format This is a data frame in the unbalanced format, that is, with one row per observation. The data consists in columns for patient identification, time of measurements, longitudinal multiple longitudinal measurements, baseline covariates, and survival data. The column names are identified as follows:
#' \describe{
#'   \item{num}{number for patient identification.}
#'   \item{sex}{gender of patient (\code{0=}Male and \code{1=}Female).}
#'   \item{age}{age of patient at day of surgery (years).}
#'   \item{time}{observed time point, with surgery date as the time origin (years).}
#'   \item{fuyrs}{maximum follow up time, with surgery date as the time origin (years).}
#'   \item{status}{censoring indicator (\code{1=}died and \code{0=}lost at follow up).}
#'   \item{grad}{valve gradient at follow-up visit.}
#'   \item{log.grad}{natural log transformation of \code{grad}.}
#'   \item{lvmi}{left ventricular mass index (standardised) at follow-up visit.}
#'   \item{log.lvmi}{natural log transformation of \code{lvmi}.}
#'   \item{ef}{ejection fraction at follow-up visit.}
#'   \item{bsa}{preoperative body surface area.}
#'   \item{lvh}{preoperative left ventricular hypertrophy.}
#'   \item{prenyha}{preoperative New York Heart Association (NYHA) classification (\code{1=}I/II and \code{3=}III/IV).}
#'   \item{redo}{previous cardiac surgery.}
#'   \item{size}{size of the valve (millimeters).}
#'   \item{con.cabg}{concomitant coronary artery bypass graft.}
#'   \item{creat}{preoperative serum creatinine (\eqn{\mu}mol/mL).}
#'   \item{dm}{preoperative diabetes.}
#'   \item{acei}{preoperative use of ace inhibitor.}
#'   \item{lv}{preoperative left ventricular ejection fraction (LVEF) (\code{1=}good, \code{2=}moderate, and \code{3=}poor).}
#'   \item{emergenc}{operative urgency (\code{0=}elective, \code{1=}urgent, and \code{3=}emergency).}
#'   \item{hc}{preoperative high cholesterol (\code{0=}absent, \code{1=}present treated, and \code{2=}present untreated).}
#'   \item{sten.reg.mix}{aortic valve haemodynamics (\code{1=}stenosis, \code{2=}regurgitation, \code{3=}mixed).}
#'   \item{hs}{implanated aortic prosthesis type (\code{1=}homograft and \code{0=}stentless porcine tissue).}
#' }
#' @source Mr Eric Lim (\url{http://www.drericlim.com})
#' @references Lim E, Ali A, Theodorou P, Sousa I, Ashrafian H, Chamageorgakis T, Duncan M, Diggle P and Pepper J. A longitudinal study of the profile and predictors of left ventricular mass regression after stentless aortic valve replacement. \emph{Ann Thorac Surg.} 2008; \strong{85(6)}: 2026-2029.
"heart.valve"
