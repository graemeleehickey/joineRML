#' Aortic valve replacement surgery data
#'
#' @description This is longitudinal data on an observational study on detecting
#'   effects of different heart valves, differing on type of tissue, implanted
#'   in the aortic position. The data consists of longitudinal measurements on
#'   three different heart function outcomes, after surgery occurred. There are
#'   several baseline covariates available, and also survival data.
#'
#' @usage data(heart.valve)
#' @format This is a data frame in the unbalanced format, that is, with one row
#'   per observation. The data consists in columns for patient identification,
#'   time of measurements, longitudinal multiple longitudinal measurements,
#'   baseline covariates, and survival data. The column names are identified as
#'   follows: \describe{
#'
#'   \item{\code{num}}{number for patient identification.}
#'
#'   \item{\code{sex}}{gender of patient (\code{0=}Male and \code{1=}Female).}
#'
#'   \item{\code{age}}{age of patient at day of surgery (years).}
#'
#'   \item{\code{time}}{observed time point, with surgery date as the time
#'   origin (years).}
#'
#'   \item{\code{fuyrs}}{maximum follow up time, with surgery date as the time
#'   origin (years).}
#'
#'   \item{\code{status}}{censoring indicator (\code{1=}died and \code{0=}lost
#'   at follow up).}
#'
#'   \item{\code{grad}}{valve gradient at follow-up visit.}
#'
#'   \item{\code{log.grad}}{natural log transformation of \code{grad}.}
#'
#'   \item{\code{lvmi}}{left ventricular mass index (standardised) at follow-up
#'   visit.}
#'
#'   \item{\code{log.lvmi}}{natural log transformation of \code{lvmi}.}
#'
#'   \item{\code{ef}}{ejection fraction at follow-up visit.}
#'
#'   \item{\code{bsa}}{preoperative body surface area.}
#'
#'   \item{\code{lvh}}{preoperative left ventricular hypertrophy.}
#'
#'   \item{\code{prenyha}}{preoperative New York Heart Association (NYHA)
#'   classification (\code{1=}I/II and \code{3=}III/IV).}
#'
#'   \item{\code{redo}}{previous cardiac surgery.}
#'
#'   \item{\code{size}}{size of the valve (millimeters).}
#'
#'   \item{\code{con.cabg}}{concomitant coronary artery bypass graft.}
#'
#'   \item{\code{creat}}{preoperative serum creatinine (\eqn{\mu}mol/mL).}
#'
#'   \item{\code{dm}}{preoperative diabetes.}
#'
#'   \item{\code{acei}}{preoperative use of ace inhibitor.}
#'
#'   \item{\code{lv}}{preoperative left ventricular ejection fraction (LVEF)
#'   (\code{1=}good, \code{2=}moderate, and \code{3=}poor).}
#'
#'   \item{\code{emergenc}}{operative urgency (\code{0=}elective, \code{1 =
#'   }urgent, and \code{3=}emergency).}
#'
#'   \item{\code{hc}}{preoperative high cholesterol (\code{0=}absent, \code{1
#'   =}present treated, and \code{2=}present untreated).}
#'
#'   \item{\code{sten.reg.mix}}{aortic valve haemodynamics (\code{1=}stenosis,
#'   \code{2=}regurgitation, \code{3=}mixed).}
#'
#'   \item{\code{hs}}{implanted aortic prosthesis type (\code{1=}homograft and
#'   \code{0=}stentless porcine tissue).}
#'
#'   }
#' @keywords datasets
#' @seealso \code{\link{pbc2}}, \code{\link{renal}},
#'   \code{\link{epileptic.qol}}.
#' @source Mr Eric Lim (\url{http://www.drericlim.com})
#' @references
#'
#' Lim E, Ali A, Theodorou P, Sousa I, Ashrafian H, Chamageorgakis T, Duncan M,
#' Diggle P, Pepper J. A longitudinal study of the profile and predictors of
#' left ventricular mass regression after stentless aortic valve replacement.
#' \emph{Ann Thorac Surg.} 2008; \strong{85(6)}: 2026-2029.
"heart.valve"
