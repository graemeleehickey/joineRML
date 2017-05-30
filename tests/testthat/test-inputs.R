library(joineRML)
context("Inputs")


test_that("too many datasets throws error", {
  # load data + fit model
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  f <- function() {
    mjoint(
      formLongFixed = list("log.bil" = log.b ~ year),
      formLongRandom = list("log.bil" = ~ 1 | id),
      formSurv = Surv(years, status2) ~ age,
      data = list(pbc2, pbc2),
      timeVar = "year",
      control = list(convCrit = "sas", rav = 0.01),
      verbose = FALSE)
  }
  # test
  expect_error(f(), "The number of datasets expected is K = 1")
})


test_that("too few datasets throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      data = list(hvd),
      inits = list("gamma" = c(0.11, 1.51, 0.80)),
      timeVar = "time",
      control = list(convCrit = "sas", rav = 0.01),
      verbose = FALSE)
  }
  # test
  expect_error(f(), "The number of datasets expected is K = 2")
})


test_that("timeVar length mismatch throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      data = hvd,
      timeVar = c("time", "time", "time"),
      control = list(convCrit = "sas", rav = 0.01),
      verbose = FALSE)
  }
  # test
  expect_error(f(), "The length of timeVar must equal 2")
})


test_that("misspelled timeVar mismatch throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      data = hvd,
      timeVar = "Time",
      control = list(convCrit = "sas", rav = 0.01),
      verbose = FALSE)
  }
  # test
  expect_error(f(), "undefined columns selected")
})


test_that("unmatched control prarameter throws warning", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      data = hvd,
      timeVar = "time",
      control = list(convCrit = "abs", tol0 = 0.1, burnin = 3,
                     mcmaxIter = 4, fake_param = 5),
      verbose = FALSE)
  }
  # test
  expect_warning(f(), "Unknown arguments passed to 'control': fake_param")
})


test_that("unmatched inits throws warning", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      inits = list("fake_param" = 5),
      data = hvd,
      timeVar = "time",
      control = list(convCrit = "abs", tol0 = 0.1, burnin = 3,
                     mcmaxIter = 4, fake_param = 5),
      verbose = FALSE)
  }
  # test
  expect_warning(f(), "Unknown initial parameters passed to 'inits': fake_param")
})


test_that("measurement time after event time throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  hvd[hvd$num == "6", "fuyrs"] <- 0.1
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      data = hvd,
      timeVar = "time",
      control = list(convCrit = "sas", rav = 0.01, burnin = 5,
                     mcmaxIter = 10),
      verbose = FALSE)
  }
  # test
  expect_error(f(), "Longitudinal measurements should not be recorded after the event time")
})


test_that("measurement time after event time throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  D <- matrix(1:9, nrow = 3)
  f <- function() {
    mjoint(
      formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                           "lvmi" = log.lvmi ~ time + sex),
      formLongRandom = list("grad" = ~ 1 | num,
                            "lvmi" = ~ time | num),
      formSurv = Surv(fuyrs, status) ~ age,
      inits = list("D" = D),
      data = hvd,
      timeVar = "time",
      control = list(convCrit = "sas", rav = 0.01, burnin = 4,
                     mcmaxIter = 5),
      verbose = FALSE)
  }
  # test
  expect_warning(f(), "Initial parameter matrix D is non positive definite: falling back to automated value")
})


test_that("unbalanced data throws warning for survival inits", {
  # load data + fit model
  data(heart.valve)
  # make unbalanced dataset for patients in common
  id <- unique(heart.valve$num[!is.na(heart.valve$log.grad)])
  hvd1 <- subset(heart.valve[!is.na(heart.valve$log.grad), ], num %in% id)
  hvd2 <- subset(heart.valve[!is.na(heart.valve$log.lvmi), ], num %in% id)
  f <- function() {mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd1, hvd2),
    timeVar = "time",
    control = list(convCrit = "sas", burnin = 5, mcmaxIter = 10),
    verbose = FALSE)
  }
  # test
  expect_message(f(), "Data are unbalanced... using sub-optimal initial parameters for gamma")
})


test_that("argument not correct object class", {
  expect_error(vcov.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(summary.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(ranef.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(plotConvergence(1), "Use only with 'mjoint' model objects.")
  expect_error(plot.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(confint.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(bootSE(1), "Use only with 'mjoint' model objects.")
  expect_error(fixef.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(formula.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(getVarCov.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(sigma.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(logLik.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(print.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(baseHaz(1), "Use only with 'mjoint' model objects.")
  expect_error(sampleData(1), "Use only with 'mjoint' model objects.")
  expect_error(print.dynSurv(1), "Use only with 'dynSurv' objects.")
  expect_error(print.dynLong(1), "Use only with 'dynLong' objects.")
  expect_error(plot.dynSurv(1), "Use only with 'dynSurv' objects.")
  expect_error(plot.dynLong(1), "Use only with 'dynLong' objects.")
  expect_error(print.summary.mjoint(1), "Use only with 'summary.mjoint' objects.")
})


test_that("argument not a summary.mjoint object", {
})


test_that("formula with unspecified longitudinal measure throws error", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = hvd,
    timeVar = "time",
    control = list(convCrit = "abs", rav = 0.05, burnin = 5,
                   mcmaxIter = 7),
    verbose = FALSE)
  # tests
  expect_error(formula(fit, process = "Longitudinal", k = NA),
               "Must specify a longitudinal outcome.")
  expect_error(formula(fit, process = "Longitudinal", k = 3),
               "Incompatible with dimensions of the joint model.")
})


test_that("simulation errors", {
  # tests
  expect_error(simData(10, model = "fake"), "Unknown model: fake")
  expect_error(simData(10, sigma2 = 1, gamma.y = 1, D = 1,
                       beta = matrix(rep(1, 4), nrow = 1)),
               "Error: this function on simulates multivariate data")
  expect_error(simData(10, D = matrix(rep(1, 16), 4, 4)),
               "Covariance matrix must be positive semi-definite")
  D <- diag(4)
  D[1, 2] <- 0.1
  expect_error(simData(10, D = D),
               "Covariance matrix is not symmetric")
})


test_that("mismatched dimensions of inits", {
  # load data
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  # tests
  expect_error(mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ year | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    inits = list("sigma2" = c(1, 2))),
    "Dimension of sigma2 inits does not match model.")
  expect_error(mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ year | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    inits = list("gamma" = c(1, 2, 3))),
    "Dimension of gamma inits does not match model.")
  expect_error(mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ year | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    inits = list("D" = diag(1, ncol = 1))),
    "Dimension of D inits does not match model.")
  expect_error(mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ year | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    inits = list("beta" = 1)),
    "Dimension of beta inits does not match model.")
})

test_that("same patients measured on all markers", {
  # load data + function to fit model
  data(heart.valve)
  hvd1 <- heart.valve[!is.na(heart.valve$log.grad), ]
  hvd2 <- heart.valve[!is.na(heart.valve$log.lvmi), ]
  fit <- function() {
    mjoint(formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                                "lvmi" = log.lvmi ~ time + sex),
           formLongRandom = list("grad" = ~ 1 | num,
                                 "lvmi" = ~ time | num),
           formSurv = Surv(fuyrs, status) ~ age,
           data = list(hvd1, hvd2),
           inits = list("gamma" = c(0.11, 1.51, 0.80)),
           timeVar = "time",
           verbose = TRUE)
  }
  # tests
  expect_error(fit(), "Every subject must have at least one measurement per each outcome")
})
