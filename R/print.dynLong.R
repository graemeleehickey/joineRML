#' @keywords internal
#' @export
print.dynLong <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "dynLong")) {
    stop("Use only with 'dynLong' objects.\n")
  }

  out <- lapply(x$pred, round, digits = digits)
  print(out)

  if (x$type == "simulated") {
    label <- if (identical(x$sampler, "conditional-rejection")) {
      "Rejection-sampling acceptance rate: "
    } else {
      "M-H acceptance rate: "
    }
    cat(paste0("\n", label, round(100 * x$accept, 1), "%"))
  }

  cat("\n")
  invisible(x)

}
