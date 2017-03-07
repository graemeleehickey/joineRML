library(joineRML)
context("Bootstrap")

test_that("bootstrap MV models", {
  skip_on_cran()
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    control = list(convCrit = "abs", tol0 = 0.05,
                   earlyPhase = 5, mcmaxIter = 100),
    verbose = FALSE)
  fit.boot <- bootSE(fit, nboot = 3)
  # tests
  expect_is(fit.boot, "bootSE")
  expect_output(str(fit.boot), "List of 11")
  expect_output(print(fit.boot))
})
