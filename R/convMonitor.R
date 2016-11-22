#' Evaluate the convergence of the MCEM algorithm
#'
#' @keywords internal
convMonitor <- function(theta, theta.new, con, verbose) {

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

  # Relative parameter change (excl. near-zero terms): for reporting only
  reldelta.pars2 <- sapply(c("beta", "sigma2", "gamma"), function(i) {
    abs(theta[[i]] - theta.new[[i]]) / (abs(theta[[i]]) + con$tol1)
  })
  max.reldelta.pars2 <- max(unlist(reldelta.pars)[mag])

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
    cat(paste("Maximum relative parameter change =",
              round(max.reldelta.pars, 6), "for",
              names(which.max(abs(unlist(reldelta.pars)))), "\n"))
    cat(paste("       --> excl. parameters <", con$rav, "=",
              round(max.reldelta.pars2, 6), "for",
              names(which.max(abs(unlist(reldelta.pars)[mag]))), "\n"))
    cat(paste("Converged:", conv, "\n\n"))
  }

  return(list(conv = conv, max.reldelta.pars = max.reldelta.pars))

}
