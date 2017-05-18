#' @keywords internal
#' @export
print.dynSurv <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "dynSurv")) {
    stop("Use only with 'dynSurv' objects.")
  }

  print(round(x$pred, digits = digits))
  cat("\n")
  invisible(x)

}
