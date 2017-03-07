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


test_that("init length mismatch throws error", {
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
      inits = list("gamma" = c(0.11, 1.51)),
      timeVar = "time",
      control = list(convCrit = "sas", rav = 0.01),
      verbose = FALSE)
  }
  # test
  expect_error(f())
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
      control = list(convCrit = "sas", rav = 0.01, earlyPhase = 5,
                     mcmaxIter = 10, fake_param = 5),
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
      control = list(convCrit = "sas", rav = 0.01, earlyPhase = 5,
                     mcmaxIter = 10),
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
      inits = list("fake_param" = 5),
      data = hvd,
      timeVar = "time",
      control = list(convCrit = "sas", rav = 0.01, earlyPhase = 5,
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
      control = list(convCrit = "sas", rav = 0.01, earlyPhase = 5,
                     mcmaxIter = 10),
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
    control = list(convCrit = "sas", earlyPhase = 5, mcmaxIter = 10),
    verbose = FALSE)
  }
  # test
  expect_message(f(), "Data are unbalanced... using sub-optimal initial parameters for gamma")
})


test_that("argument not an mjoint object", {
  expect_error(vcov.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(summary.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(ranef.mjoint(1), "Use only with 'mjoint' model objects.")
  expect_error(plotConvergence(1), "Use only with 'mjoint' model objects.")
  expect_error(plot.mjoint(1), "Use only with 'mjoint' model objects.")
})
