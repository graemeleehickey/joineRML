#' @keywords internal
convMonitor <- function(theta, theta.new, log.lik, log.lik.new, con, verbose) {

  # Absolute parameter change
  absdelta.pars <- sapply(c("D", "beta", "sigma2", "gamma"), function(i) {
    theta[[i]] - theta.new[[i]]
  })
  max.absdelta.pars <- max(abs(unlist(absdelta.pars)))
  cond1 <- (max.absdelta.pars < con$tol0)

  # Relative parameter change
  reldelta.pars <- sapply(c("D", "beta", "sigma2", "gamma"), function(i) {
    abs(theta[[i]] - theta.new[[i]]) / (abs(theta[[i]]) + con$tol1)
  })
  max.reldelta.pars <- max(unlist(reldelta.pars))
  cond2 <- (max.reldelta.pars < con$tol2)

  # Either parameter change satisfied
  cond3 <- (cond1 || cond2)

  # SAS criteria
  mag <- (abs(unlist(theta[-which(names(theta) == "haz")])) > con$rav)
  sasdelta <- function() {
    d1 <- unlist(absdelta.pars) < con$tol0 # abs
    d2 <- unlist(reldelta.pars) < con$tol2 # rel
    all(d1[!mag]) & all(d2[mag])
  }
  cond4 <- sasdelta()

  # Absolute parameter change (excl. large terms): for reporting only
  max.absdelta.pars2 <- ifelse(any(!mag), max(unlist(absdelta.pars)[!mag]), NA)

  # Relative parameter change (excl. near-zero terms): for reporting only
  max.reldelta.pars2 <- ifelse(any(mag), max(unlist(reldelta.pars)[mag]), NA)

  # Log-likelihood: for reporting only
  rel.ll <- (log.lik.new - log.lik) / abs(log.lik + con$tol1)

  # Choose convergence criterion to use
  if (con$convCrit == "abs") {
    conv <- cond1
  } else if (con$convCrit == "rel") {
    conv <- cond2
  } else if (con$convCrit == "either") {
    conv <- cond3
  } else if (con$convCrit == "sas") {
    conv <- cond4
  }

  # Monitoring messages
  if (verbose) {
    cat(paste("Maximum absolute parameter change =",
              round(max.absdelta.pars, 6), "for",
              names(which.max(abs(unlist(absdelta.pars)))), "\n"))
    if (any(!mag)) {
      cat(paste("      ---> for parameters  <", con$rav, "=",
                round(max.absdelta.pars2, 6), "for",
                names(which.max(abs(unlist(absdelta.pars)[!mag]))), "\n"))
    }
    cat(paste("Maximum relative parameter change =",
              round(max.reldelta.pars, 6), "for",
              names(which.max(abs(unlist(reldelta.pars)))), "\n"))
    if (any(mag)) {
      cat(paste("      ---> for parameters >=", con$rav, "=",
                round(max.reldelta.pars2, 6), "for",
                names(which.max(unlist(reldelta.pars)[mag])), "\n"))
    }
    cat(paste("Relative change in log-likelihood =", round(rel.ll, 6), "\n"))
    cat(paste("Converged:", conv, "\n\n"))
  }

  return(list("conv" = conv,
              "max.reldelta.pars" = max.reldelta.pars,
              "rel.ll" = rel.ll))

}
