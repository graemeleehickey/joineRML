#' @keywords internal
#' @export
print.mjoint <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  K <- x$dims$K

  #*****************************************************
  # Call statements
  #*****************************************************

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  #*****************************************************
  # Data descriptives
  #*****************************************************

  cat("Number of subjects:", x$dims$n, "\n")
  cat("Number of longitudinal outcomes: K =", x$dims$K, "\n")
  cat("Number of observations:\n")
  for (k in 1:K) {
    if(any(is.null(names(x$formLongFixed)))) {
      cat("   Outcome ", k, ": n = ", x$dims$nk[k], "\n", sep = "")
    } else {
      cat("   Outcome ", k, " (", names(x$formLongFixed)[k], "): n = ", x$dims$nk[k],
          "\n", sep = "")
    }
  }

  #*****************************************************
  # Model summaries
  #*****************************************************

  cat("\nJoint Model Summary:\n")
  if (K == 1) {
    cat("\nLongitudinal Process: Univariate linear mixed-effects model\n")
  } else {
    cat("\nLongitudinal Process: Multivariate linear mixed-effects model\n")
  }
  for (k in 1:K) {
    cat("    ", paste0(deparse(x$formLongFixed[[k]]), ", random = ",
                       deparse(x$formLongRandom[[k]]), "\n"))
  }

  cat("Event Process: Cox proportional hazards model\n")
  cat("     ", paste(deparse(x$formSurv), sep = "\n", collapse = "\n"), "\n", sep = "")

  #*****************************************************
  # Variance components
  #*****************************************************

  cat("\nVariance Components:\n\n")

  # Random effects variance-covariance matrix
  print(getVarCov(x))

  #*****************************************************
  # Model coefficients
  #*****************************************************

  cat("\nCoefficient Estimates:\n")

  print(coef(x)[c("beta", "gamma")])

  cat("\n")
  invisible(x)

}
