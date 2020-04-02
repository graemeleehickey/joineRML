#' Sample from an \code{mjoint} object
#'
#' @description Generic function used to sample a subset of data from an object
#'   of class \code{mjoint} with a specific number of subjects.
#'
#' @inheritParams confint.mjoint
#' @param size number of subjects to include in the sampled subset. If
#'   \code{size=NULL} (default), then size is set equal to the number of
#'   subjects used to fit the \code{mjoint} model.
#' @param replace use replacement when sampling subjects? Default is
#'   \code{TRUE}. If replacement is used, then the subjects are re-labelled from
#'   1 to \code{size}.
#'
#' @details This function is primarily intended for internal use in the
#'   \code{\link{bootSE}} function in order to permit bootstrapping. However, it
#'   can be used for other purposes given a fitted \code{mjoint} object.
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
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
sampleData <- function(object, size = NULL, replace = TRUE) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  # Extract from model fit
  data <- object$data
  survData <- object$survData
  id <- object$id
  id.labs <- unique(survData[ , id])
  surv.time.lab <- all.vars(object$formSurv)[1]
  n <- length(id.labs)
  K <- object$dims$K
  size <- ifelse(is.null(size), n, size)

  if (!replace) {
    if (!is.null(size)) {
      if (size > n) {
        stop("Cannot select more subjects than in data without replacement")
      }
    }
  }

  # Random sample of subjects (with replacement)
  i <- sample(id.labs, size = size, replace = replace)

  # Longitudinal data
  longData.boot <- list()
  for (k in 1:K) {
    out <- lapply(i, function(u) data[[k]][data[[k]][, id] == u, ])
    m <- do.call("c", lapply(out, nrow))
    longData.boot[[k]] <- do.call("rbind", out)
    if (replace) {
      id.new <- rep(1:size, m)
      longData.boot[[k]][, id] <- id.new
    }
  }

  # Time-to-event data
  survData.boot <- survData[match(i, survData[ , id]), ]
  if (replace) {
    survData.boot[, id] <- 1:size
  }
  survData.boot[, surv.time.lab] <- survData.boot[, surv.time.lab]

  return(list(longData.boot = longData.boot,
              survData.boot = survData.boot))

}
