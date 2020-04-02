#' @keywords internal
process_newdata <- function(object, newdata, newSurvData = NULL, tobs = NULL) {

  if (!inherits(object, "mjoint")) {
    stop("Use only with 'mjoint' model objects.\n")
  }

  K <- object$dims$K
  r <- object$dims$r
  p <- object$dims$p
  q <- object$dims$q

  if (class(newdata) != "list") {
    balanced <- TRUE
    newdata <- list(newdata)
    if (K > 1) {
      for (k in 2:K) {
        newdata[[k]] <- newdata[[1]]
      }
    }
  } else {
    balanced <- (length(unique(newdata)) == 1)
  }
  if (length(newdata) != K) {
    stop(paste("The number of datasets expected is K =", K))
  }

  #*****************************************************
  # Multivariate longitudinal data
  #*****************************************************

  Xk.new <- list()
  tk.new <- list()
  Zk.new <- list()
  yk.new <- list()
  ffk <- list()
  nk <- vector(length = K)

  for (k in 1:K) {
    termsX <- object$lfit[[k]]$terms
    origData <- model.frame(termsX, object$data[[k]])
    xlev <- .getXlevels(termsX, origData)
    mfX.new <- model.frame(termsX, newdata[[k]], xlev = xlev)
    ffk[[k]] <- nlme::splitFormula(object$formLongRandom[[k]], "|")[[1]]

    Xk.new[[k]] <- model.matrix(termsX, mfX.new)
    tk.new[[k]] <- newdata[[k]][[object$timeVar[[k]]]]
    yk.new[[k]] <- model.response(mfX.new, "numeric")
    Zk.new[[k]] <- model.matrix(ffk[[k]], newdata[[k]])
    nk[k] <- nrow(Xk.new[[k]])
  }

  X.new <- as.matrix(Matrix::bdiag(Xk.new))
  Z.new <- as.matrix(Matrix::bdiag(Zk.new))
  y.new <- unlist(yk.new)

  #*****************************************************
  # Time-to-event data
  #*****************************************************

  if (q > 0) {
    termsT <- object$sfit$terms
    if (is.null(newSurvData)) {
      mfT.new <- model.frame(delete.response(termsT), newdata[[1]],
                             xlev = object$sfit$xlevels)
    } else {
      mfT.new <- model.frame(delete.response(termsT), newSurvData,
                             xlev = object$sfit$xlevels)
    }
    v.new <- model.matrix(delete.response(termsT), mfT.new)[1, -1]
    v.new <- v.new - object$dmats$t$xcenter
  } else {
    v.new <- NULL
  }

  survdat2 <- object$dmats$t$survdat2[object$dmats$t$survdat2$delta == 1, ]
  survdat2 <- survdat2[order(survdat2$T), ]
  tmax <- max(survdat2$T)

  if (is.null(tobs)) {
    tobs <- 0
    for (k in 1:K) {
      tobs <- max(tobs, max(newdata[[k]][, object$timeVar[k]]))
    }
  } else {
    if (tobs > tmax) {
      stop("Cannot extrapolate beyond final failure time")
    }
  }

  tj <- object$dmats$t$tj # if tj.ind = 0, deal with it elsewhere
  tj.ind <- sum(tj <= tobs)

  #*****************************************************
  # Z(t_j) for unique failure times t_j
  #*****************************************************

  Zdat.fail <- data.frame(tj[1:max(tj.ind, 1)])
  Zk.fail <- list()

  for (k in 1:K) {
    names(Zdat.fail) <- object$timeVar[[k]]
    Zk.fail[[k]] <- model.matrix(ffk[[k]], Zdat.fail)
  }

  Z.fail <- as.matrix(Matrix::bdiag(Zk.fail))
  IW.fail <- do.call("cbind", lapply(1:K, function(i) diag(max(tj.ind, 1))))

  #*****************************************************
  # Output
  #*****************************************************

  out <- list(
    X = X.new,
    Xk = Xk.new,
    tk = tk.new,
    Z = Z.new,
    Zk = Zk.new,
    y = y.new,
    yk = yk.new,
    nk = nk,
    v = v.new,
    tobs = tobs,
    tmax = tmax,
    tj.ind = tj.ind,
    Z.fail = Z.fail,
    IW.fail = IW.fail,
    K = K, p = p, r = r, q = q
  )

  return(out)

}
