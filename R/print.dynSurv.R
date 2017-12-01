#' @keywords internal
#' @export
print.dynSurv <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "dynSurv")) {
    stop("Use only with 'dynSurv' objects.\n")
  }

  if (!is.null(x$horizon)) {
    cat(paste0("Last follow-up time = ", round(x$data.t$tobs, digits), "\n"))
    cat(paste0("       Horizon time = ", round(x$horizon, digits), "\n\n"))
  }

  print(round(x$pred, digits = digits))

  if (x$type == "simulated") {
    cat(paste0("\nM-H acceptance rate: ", round(100 * x$accept, 1), "%"))
  }

  cat("\n")
  invisible(x)

}
