% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getVarCov.mjoint.R
\name{getVarCov.mjoint}
\alias{getVarCov.mjoint}
\title{Extract variance-covariance matrix of random effects from an \code{mjoint}
object}
\usage{
\method{getVarCov}{mjoint}(obj, ...)
}
\arguments{
\item{obj}{an object inheriting from class \code{mjoint} for a joint model of
time-to-event and multivariate longitudinal data.}

\item{...}{additional arguments; currently none are used.}
}
\value{
A variance-covariance matrix.
}
\description{
Extract variance-covariance matrix of random effects from an
  \code{mjoint} object.
}
\examples{
\dontrun{
# Fit a joint model with bivariate longitudinal outcomes

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

fit2 <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    verbose = TRUE)

getVarCov(fit2)
}
}
\references{
Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
Springer Verlag; 2000.
}
\seealso{
\code{\link[nlme]{getVarCov}} for the generic method description.
}
\author{
Graeme L. Hickey (\email{graemeleehickey@gmail.com})
}
\keyword{methods}
