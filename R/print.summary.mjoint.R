#' @keywords internal
#' @export
print.summary.mjoint <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "summary.mjoint")) {
    stop("Use only with 'summary.mjoint' objects.\n")
  }

  K <- x$dims$K

  #*****************************************************
  # Call statements
  #*****************************************************

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  #*****************************************************
  # Data descriptives
  #*****************************************************

  event.pc <- with(x$sfit, round(100 * nevent / n, 1))

  cat("\nData Descriptives:\n")

  cat("\nEvent Process\n")
  cat("    Number of subjects:", x$dims$n, "\n")
  cat("    Number of events: ", x$sfit$nevent, " (", event.pc, "%)\n", sep = "")

  cat("\nLongitudinal Process\n")
  cat("    Number of longitudinal outcomes: K =", x$dims$K, "\n")
  cat("    Number of observations:\n")
  for (k in 1:K) {
    if(any(is.null(names(x$formLongFixed)))) {
      cat("      Outcome ", k, ": n = ", x$dims$nk[k], "\n", sep = "")
    } else {
      cat("      Outcome ", k, " (", names(x$formLongFixed)[k], "): n = ", x$dims$nk[k],
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

  cat("Model fit statistics:\n")
  model.sum <- data.frame(log.Lik = x$logLik, AIC = x$AIC, BIC = x$BIC, row.names = "")
  print(model.sum)

  #*****************************************************
  # Variance components
  #*****************************************************

  cat("\nVariance Components:\n\n")

  # Random effects variance-covariance matrix
  print(x$D)

  sigma <- x$sigma
  cat("\nResidual standard errors(s):\n")
  print(sigma)

  #*****************************************************
  # Model coefficients
  #*****************************************************

  cat("\nCoefficient Estimates:\n")

  cat("\nLongitudinal sub-model:\n")
  out <- as.data.frame(round(x$"coefs.long", digits))
  ind <- out$"p-value" == 0
  out$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), out$"p-value")
  out$"p-value"[ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
  print(out)

  cat("\nTime-to-event sub-model:\n")
  out <- as.data.frame(round(x$"coefs.surv", digits))
  ind <- out$"p-value" == 0
  out$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), out$"p-value")
  out$"p-value"[ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
  print(out)

  #*****************************************************
  # Computational statistics
  #*****************************************************

  cat("\nAlgorithm Summary:\n")

  cat("    Total computational time:", round(x$comp.time[1], 1),
      attr(x$comp.time[1], "units"), "\n")
  cat("    EM algorithm computational time:", round(x$comp.time[2], 1),
      attr(x$comp.time[2], "units"), "\n")
  cat("    Convergence status:", ifelse(x$conv, "converged\n", "failed\n"))
  cat("    Convergence criterion:", x$control$convCrit, "\n")
  cat("    Final Monte Carlo sample size:", x$finalnMC, "\n")
  cat("    Standard errors calculated using method:", x$se.type)
  if (x$se.type == "boot") {
    cat("\n    Number of bootstraps: B =", x$nboot, "\n")
    cat("    Bootstrap computational time:", round(x$boot.time, 1),
        attr(x$boot.time, "units"))
  }

  cat("\n")
  invisible(x)

}
