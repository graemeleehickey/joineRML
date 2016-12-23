#' Compute the nearest positive definite matrix to an approximate one
#'
#' @keywords internal
fast_nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                    maxits = 100) {
  # Code is from JM package (v. 1.4-5) by Dimitris Rizopoulos
  # Copied here due to function being non-exported by JM package, and also because
  # the Matrix package version is an order of magnitude slower.
  # This code is originally based on function nearcor() submitted to R-help by
  # Jens Oehlschlagel on 2007-07-13, and function posdefify() from package 'sfsmisc'
  if (!(is.numeric(M) && is.matrix(M))) {
    stop("Input matrix M must be square and symmetric.\n")
  }
  inorm <- function (x) {
    max(rowSums(abs(x)))
  }
  n <- ncol(M)
  U <- matrix(0, n, n)
  X <- M
  iter <- 0
  converged <- FALSE
  while (iter < maxits && !converged) {
    Y <- X
    T <- Y - U
    e <- eigen(Y, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    D <- if (length(d) > 1) {
      diag(d)
    } else {
      as.matrix(d)
    }
    p <- (d > eig.tol * d[1])
    QQ <- Q[, p, drop = FALSE]
    X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
    U <- X - T
    X <- (X + t(X)) / 2
    conv <- inorm(Y - X) / inorm(Y)
    iter <- iter + 1
    converged <- conv <= conv.tol
  }
  X <- (X + t(X)) / 2
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag) / diag(X))
    X[] <- D * X * rep(D, each = n)
  }
  (X + t(X)) / 2
}
