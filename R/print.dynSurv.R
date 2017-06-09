#' @keywords internal
#' @export
print.dynSurv <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "dynSurv")) {
    stop("Use only with 'dynSurv' objects.\n")
  }

  print(round(x$pred, digits = digits))

  if (x$type == "simulated") {
    cat(paste0("\nM-H acceptance rate: ", round(100 * x$accept, 1), "%"))
  }

  cat("\n")
  invisible(x)

}
