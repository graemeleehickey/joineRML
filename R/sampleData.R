#' Sample from an \code{mjoint} object
#'
#' Generic function used to sample a subset of data from an object of class
#' \code{mjoint}, with a specific size of number of subjects.
#'
#' @inheritParams confint.mjoint
#' @param size number of subjects to include in the sampled subset. If
#'   \code{size=NULL} (default), then size is set equal to the number of
#'   subjects used to fit the \code{mjoint} model.
#'
#' @details This function is primarily intended for internal use in the
#'   \code{\link{bootSE}} function in order to permit bootstrapping. However, it
#'   can be used for other purposes given a fitted \code{mjoint} object.
#'
#' @author Graeme L. Hickey (\email{graeme.hickey@@liverpool.ac.uk})
#' @keywords methods datagen multivariate survival
#' @seealso \code{\link{mjoint}}.
#'
#' @return A list of 2 data.frames: one recording the requisite longitudinal
#'   outcomes data, and the other recording the time-to-event data.
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
#' sampleData(fit2, size = 10)
#' }
sampleData <- function(object, size = NULL) {

  # Extract from model fit
  data <- object$data
  survData <- object$survData
  id <- object$id
  id.labs <- unique(survData[ , id])
  surv.time.lab <- all.vars(object$formSurv)[1]
  n <- length(id.labs)
  K <- object$dims$K

  # Random sample of subjects (with replacement)
  i <- sample(id.labs, size = ifelse(is.null(size), n, size),
              replace = TRUE)

  # Longitudinal data
  longData.boot <- list()
  for (k in 1:K) {
    out <- lapply(i, function(u) data[[k]][data[[k]][, id] == u, ])
    m <- do.call("c", lapply(out, nrow))
    id.new <- rep(1:n, m)
    longData.boot[[k]] <- do.call("rbind", out)
    longData.boot[[k]][, id] <- id.new
  }

  # Time-to-event data
  survData.boot <- survData[match(i, survData[ , id]), ]
  survData.boot[ , id] <- 1:n
  survData.boot[surv.time.lab] <- survData.boot[surv.time.lab]

  return(list(longData.boot = longData.boot,
              survData.boot = survData.boot))

}
