library(joineRML)
context("Models fit")


test_that("univariate random-intercept model works", {
  # load data + fit model
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  fit <- mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "sas", rav = 0.01),
    verbose = FALSE)
  # tests
  expect_is(fit, "mjoint")
  expect_true(fit$conv)
  expect_equal(length(fixef(fit)), 2)
  expect_equal(nrow(ranef(fit)), fit$dims$n)
  expect_output(str(coef(fit)), "List of 5")
})

test_that("multivariate model works", {
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
    control = list(convCrit = "sas", rav = 0.01),
    verbose = TRUE)
  fit.summ <- summary(fit)
  # tests
  expect_is(fit, "mjoint")
  expect_true(fit$conv)
  expect_equal(length(fixef(fit)), 7)
  expect_equal(length(fixef(fit, process = "Event")), 3)
  expect_equal(nrow(ranef(fit)), fit$dims$n)
  expect_equal(dim(attr(ranef(fit, postVar = TRUE), "postVar")), c(3, 3, fit$dims$n))
  expect_output(str(coef(fit)), "List of 5")
  expect_equal(coef(fit)$gamma, c(0.1088112, 1.5183319, 0.7971334),
               tolerance = 0.05, check.attributes = FALSE)
  expect_output(print(fit))
  expect_output(str(sigma(fit)), "Named num")
  expect_output(str(confint(fit)), "num")
  expect_equal(dim(confint(fit)), c(10, 2))
  expect_is(getVarCov(fit), c("random.effects", "VarCov"))
  expect_output(str(AIC(fit)), "num")
  expect_output(str(logLik(fit)), "Class 'logLik'")
  expect_is(fit.summ, "summary.mjoint")
  expect_output(str(fit.summ), "List of 20")
  expect_output(print(fit.summ))
  expect_is(vcov(fit), "matrix")
  expect_equal(dim(vcov(fit)), c(18, 18))
  expect_is(vcov(fit, correlation = TRUE), "matrix")
  expect_equal(dim(vcov(fit, correlation = TRUE)), c(18, 18))
  expect_is(formula(fit, process = "Longitudinal"), "formula")
  expect_is(formula(fit, process = "Event"), "formula")
  expect_equal(formula(fit, process = "Longitudinal"),
               formula(log.grad ~ time + sex + hs))
})


test_that("different convergence criteria work", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  fit1 <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    control = list(convCrit = "sas", earlyPhase = 5, mcmaxIter = 10),
    verbose = FALSE)
  fit2 <- update(fit1, control = list(convCrit = "rel", earlyPhase = 5,
                                      mcmaxIter = 10))
  fit3 <- update(fit1, control = list(convCrit = "abs", earlyPhase = 5,
                                      mcmaxIter = 10))
  fit4 <- update(fit1, control = list(convCrit = "either", earlyPhase = 5,
                                      mcmaxIter = 10))
  # tests
  # SAS
  expect_is(fit1, "mjoint")
  expect_equal(fit1$control$convCrit, "sas")
  # rel
  expect_is(fit2, "mjoint")
  expect_equal(fit2$control$convCrit, "rel")
  # abs
  expect_is(fit3, "mjoint")
  expect_equal(fit3$control$convCrit, "abs")
  # either
  expect_is(fit4, "mjoint")
  expect_equal(fit4$control$convCrit, "either")
})


test_that("models fit to unbalanced data", {
  # load data + fit model
  data(heart.valve)
  # make unbalanced dataset for patients in common
  id <- unique(heart.valve$num[!is.na(heart.valve$log.grad)])
  hvd1 <- subset(heart.valve[!is.na(heart.valve$log.grad), ], num %in% id)
  hvd2 <- subset(heart.valve[!is.na(heart.valve$log.lvmi), ], num %in% id)
  fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd1, hvd2),
    timeVar = "time",
    control = list(convCrit = "sas", earlyPhase = 5, mcmaxIter = 10),
    verbose = FALSE)
  expect_false(nrow(hvd1) == nrow(hvd2))
  expect_is(fit, "mjoint")
})