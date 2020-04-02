#' Plot a \code{ranef.mjoint} object
#'
#' @description Displays a plot of the BLUPs and approximate 95\% prediction
#'   interval for each subject.
#'
#' @param x an object inheriting from class \code{ranef.mjoint}, representing
#'   the estimated random effects for the \code{mjoint} object from which it was
#'   produced.
#' @inheritParams confint.mjoint
#'
#' @author Graeme L. Hickey (\email{graemeleehickey@@gmail.com})
#' @keywords methods
#' @seealso \code{\link{ranef.mjoint}}.
#'
#' @references
#'
#' Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
#' Springer Verlag; 2000.
#'
#' @return an object inheriting from class \code{ggplot}, which displays a
#'   trellis plot with a separate panel for each effect, showing a dotplot (with
#'   optional error bars indicating approximate 95\% prediction intervals if the
#'   argument \code{postVar=TRUE} is set in the call to \code{\link{ranef}}) for
#'   each subject (by row).
#' @import ggplot2
#' @importFrom utils stack
#' @export
#'
#' @examples
#' \dontrun{
#' require(ggplot2)
#' data(heart.valve)
#' hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
#' set.seed(1)
#'
#' fit1 <- mjoint(formLongFixed = log.lvmi ~ time,
#'     formLongRandom = ~ time | num,
#'     formSurv = Surv(fuyrs, status) ~ 1,
#'     data = hvd,
#'     timeVar = "time")
#'
#' plot(ranef(fit1, postVar = TRUE))
#' }
plot.ranef.mjoint <- function(x, ...) {

  if (!inherits(x, "ranef.mjoint")) {
    stop("Use only with 'ranef.mjoint' objects.\n")
  }

  xstk <- utils::stack(x)
  xstk$subject <- rep(rownames(x), ncol(x))

  ranef.vars <- attr(x, "postVar")
  if (!is.null(ranef.vars)) {
    if (is.array(ranef.vars)) {
      ses <- sqrt(apply(ranef.vars, 3, diag))
      ses <- as.data.frame(t(ses))
      ses <- utils::stack(ses)
      xstk$se <- ses$values
    } else {
      xstk$se <- sqrt(ranef.vars)
    }
    xstk$xmin <- with(xstk, values - 1.96 * se)
    xstk$xmax <- with(xstk, values + 1.96 * se)
  }

  p <- ggplot(aes_string(x = 'values', y = 'subject'), data = xstk) +
    geom_point() +
    facet_grid(~ ind, scales = "free_x")

  if (!is.null(xstk$se)) {
    p <- p + geom_errorbarh(aes_string(xmin = 'xmin', xmax = 'xmax'),
                            height = 0)
  }

  return(p)

}
