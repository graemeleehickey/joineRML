#' @keywords internal
#' @export
print.dynLong <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "dynLong")) {
    stop("Use only with 'dynLong' objects.\n")
  }

  out <- lapply(x$pred, round, digits = digits)
  print(out)
  cat("\n")
  invisible(x)

}
