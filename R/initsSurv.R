#' @keywords internal
initsSurv_balanced <- function(data, lfit, survdat2, formSurv, id, timeVar, K, q) {

  W <- matrix(nrow = nrow(data[[1]]), ncol = K)
  for (k in 1:K) {
    XbetaZb <- predict(lfit[[k]], level = 1)
    Xbeta <- predict(lfit[[k]], level = 0)
    W[, k] <- XbetaZb - Xbeta
  }
  colnames(W) <- paste0("gamma_", 1:K)

  dataAG <- by(data[[1]], data[[1]][ , id], FUN = function(u) {
    id.col <- u[ , id]
    T <- survdat2[survdat2$id == id.col[1], "T"] + 1e-06
    start <- u[ , timeVar[k]]
    stop <- c(u[, timeVar[k]][-1], T)
    status <- rep(0, length(id.col))
    status[length(id.col)] <- survdat2[survdat2$id == id.col[1], "delta"]
    if (q > 0) {
      X <- survdat2[survdat2$id == id.col[1], 2:(q + 1), drop = FALSE]
      X <- X[rep(1, length(id.col)), ]
      if (q == 1) {
        X <- matrix(X, ncol = 1)
      }
      colnames(X) <- names(survdat2)[2:(q + 1)]
    }
    if (q > 0) {
      data.frame("id" = id.col, start, stop, status, X)
    } else {
      data.frame("id" = id.col, start, stop, status)
    }
  })
  dataAG <- do.call("rbind", dataAG)
  dataAG <- cbind(dataAG, W)

  formK <- paste0("gamma_", 1:K, collapse = " + ")
  if (q > 0) {
    formX <- paste0(names(survdat2)[2:(q + 1)], collapse = " + ")
    formS <- paste("Surv(start, stop, status) ~ ", formX, "+", formK)
  } else {
    formS <- paste("Surv(start, stop, status) ~", formK)
  }
  fitAG <- survival::coxph(as.formula(formS), data = dataAG)

  gamma <- coef(fitAG)
  haz <- survival::coxph.detail(fitAG)$hazard

  return(list("gamma" = gamma, "haz" = haz))

}


#' @keywords internal
initsSurv_unbalanced <- function(sfit, K, q) {

  # Internal function for generating initial parameters for the survival
  # sub-model when the data are *unbalanced*, i.e. longitudinal outcomes can be
  # measured at different times and frequencies within a patient

  # hazard
  if (q > 0) {
    haz <- survival::coxph.detail(sfit)$hazard
  } else {
    sfit.start <- survival::survfit(sfit)
    haz <- with(sfit.start, n.event / n.risk)
    haz <- haz[sfit.start$n.event > 0]
  }

  # gamma
  if (q > 0) {
    gamma.v <- coef(sfit)
  } else {
    gamma.v <- NULL
  }
  gamma.y <- rep(0, K)
  names(gamma.y) <- paste0(rep("gamma", K), "_", 1:K)
  gamma <- c(gamma.v, gamma.y)

  return(list("gamma" = gamma, "haz" = haz))

}


#' @keywords internal
initsSurv <- function(data, lfit, sfit, survdat2, formSurv, id, timeVar, K, q,
                      balanced, inits) {

  # Internal function for generating initial parameters for the survival
  # sub-model when the data are *unbalanced*, i.e. longitudinal outcomes can be
  # measured at different times and frequencies within a patient

  if (balanced & !("gamma" %in% names(inits))) {
    inits.surv <- initsSurv_balanced(
      data = data, lfit = lfit, survdat2 = survdat2, formSurv = formSurv,
      id = id, timeVar = timeVar, K = K, q = q)
  } else {
    if (!("gamma" %in% names(inits))) {
      message("Data are unbalanced... using sub-optimal initial parameters for gamma")
    }
    inits.surv <- initsSurv_unbalanced(sfit = sfit, K = K, q = q)
  }

  # over-ride with user-specified inits
  if ("gamma" %in% names(inits)) {
    gamma <- inits$gamma
    if (length(gamma) != (q + K)) {
      stop("Dimension of gamma inits does not match model.")
    }
    names(gamma) <- names(inits.surv[["gamma"]])
    inits.surv[["gamma"]] <- gamma
  }
  if ("haz" %in% names(inits)) {
    inits.surv[["haz"]] <- inits$haz
  }

  return(inits.surv)

}
