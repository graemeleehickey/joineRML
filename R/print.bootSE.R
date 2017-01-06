#' @keywords internal
#' @export
print.bootSE <- function(x, digits = max(4, getOption("digits") - 4), ...) {

  if (!inherits(x, "bootSE")) {
    stop("Use only with 'bootSE' objects.\n")
  }

  ci <- as.numeric(100*x$ci)
  if ((ci %% 1) == 0) {
    ci <- as.integer(ci)
  }

  cat("\nBootstrap SE estimates and percentile (", ci, "%) confidence intervals\n\n",
      sep = "")

  coefs.beta <- cbind("Value" = x$coefficients$beta,
                      "Std.Err" = x$beta.se,
                      "CI.lower" = x$beta.ci[1, ],
                      "CI.upper" = x$beta.ci[2, ])

  coefs.gamma <- cbind("Value" = x$coefficients$gamma,
                       "Std.Err" = x$gamma.se,
                       "CI.lower" = x$gamma.ci[1, ],
                       "CI.upper" = x$gamma.ci[2, ])

  coefs.sigma2 <- cbind("Value" = x$coefficients$sigma2,
                        "Std.Err" = x$sigma2.se,
                        "CI.lower" = x$sigma2.ci[1, ],
                        "CI.upper" = x$sigma2.ci[2, ])

  out <- rbind(coefs.beta, coefs.gamma, coefs.sigma2)
  out <- round(out, digits)
  out <- as.data.frame(out)
  pars <- rownames(out)
  submodel <- rep("", nrow(out))
  submodel[1] <- "Longitudinal"
  submodel[nrow(out) - nrow(coefs.sigma2) + 1] <- "Residual SE"
  submodel[nrow(coefs.beta) + 1] <- "Time-to-event"
  out <- cbind(submodel, pars, out)
  colnames(out) <- c("", "Coefficient", "Estimate", "SE",
                     paste0(ci, "% CI Lower"), paste0(ci, "% CI Upper"))
  rownames(out) <- NULL

  print(out, row.names = FALSE)
  cat("\nBootstrap computational time:", round(x$boot.time, 1),
      attr(x$boot.time, "units"))
  cat("\nBootstrap model convergence rate: ",
      round(100 * x$conv / x$nboot, 1), "%", sep = "")
  cat("\n")

  invisible(x)

}
