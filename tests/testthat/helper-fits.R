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
