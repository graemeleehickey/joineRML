library(joineRML)
context("Bootstrap")

test_that("bootstrap MV models", {
  skip_on_cran()
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  set.seed(123)
  fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    control = list(convCrit = "abs", tol0 = 0.1, tol.em = 1e-02,
                   earlyPhase = 5, mcmaxIter = 20),
    verbose = FALSE)
  fit.boot <- bootSE(fit, nboot = 2, verbose = TRUE)
  # tests
  expect_is(fit.boot, "bootSE")
  expect_output(str(fit.boot), "List of 11")
  expect_output(print(fit.boot))
  expect_is(summary(fit, bootSE = fit.boot), "summary.mjoint")
  expect_output(str(summary(fit, bootSE = fit.boot)), "List of 22")
  expect_equal(summary(fit, bootSE = fit.boot)$se.type, "boot")
  expect_equal(dim(confint(fit, bootSE = fit.boot)), c(10, 2))
  expect_equal(dim(confint(fit, parm = "Longitudinal", bootSE = fit.boot)), c(7, 2))
  expect_equal(dim(confint(fit, parm = "Event", bootSE = fit.boot)), c(3, 2))
})
