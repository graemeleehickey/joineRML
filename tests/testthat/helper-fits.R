# Shared model-fitting helpers, sourced automatically by testthat before tests.

# A small, fast-converging univariate joint model on the pbc2 data, used across
# the dynamic-prediction, plotting, and accessor-method tests.
fit_pbc_uni <- function() {
  data(pbc2)
  set.seed(4821)
  mjoint(
    formLongFixed = list("serBilir" = log(serBilir) ~ year),
    formLongRandom = list("serBilir" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 0.5, burnin = 20, mcmaxIter = 100),
    verbose = FALSE)
}

# Single-subject newdata for dynamic-prediction tests.
pbc_newdata <- function(subject = 1) {
  data(pbc2)
  droplevels(subset(pbc2, id == subject))
}

# Small data with equal visit counts for checking dimensions and fitted quantities.
small_joint_data <- function() {
  set.seed(73)
  n <- 30
  b1 <- rnorm(n)
  b2 <- 0.5 * b1 + rnorm(n, sd = 0.7)
  data.frame(
    id = rep(seq_len(n), each = 3),
    time = rep(0:2, n),
    y1 = 2 + rep(0:2, n) * 0.2 + rep(b1, each = 3) + rnorm(3 * n, sd = 0.4),
    y2 = 1 - rep(0:2, n) * 0.1 + rep(b2, each = 3) + rnorm(3 * n, sd = 0.4),
    T = rep(3 + rexp(n), each = 3),
    status = rep(rep(c(1, 1, 0), length.out = n), each = 3)
  )
}

fit_small_joint <- function(data = small_joint_data(),
                            formLongFixed = list(y1 = y1 ~ time),
                            formLongRandom = list(y1 = ~ 1 | id),
                            pfs = TRUE, mcmaxIter = 4,
                            formSurv = Surv(T, status) ~ 1) {
  # A fixed iteration count and fixed integration nodes make these checks fast
  # and reproducible. These settings are not intended for statistical analysis.
  survdat <- if (is.data.frame(data)) data else data[[1]]
  survmat <- model.matrix(formSurv, survdat)
  q <- sum(attr(survmat, "assign") != 0)
  suppressMessages(mjoint(
    formLongFixed = formLongFixed,
    formLongRandom = formLongRandom,
    formSurv = formSurv,
    data = data,
    timeVar = "time",
    inits = list(gamma = rep(0.1, length(formLongFixed) + q)),
    pfs = pfs,
    control = list(type = "halton", nMC = 128, nMCmax = 128,
                   burnin = 1, mcmaxIter = mcmaxIter, convCrit = "abs",
                   tol0 = Inf, tol.em = 1e-3)
  ))
}
