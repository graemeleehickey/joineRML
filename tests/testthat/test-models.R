library(joineRML)
context("Models fit")


test_that("univariate random-intercept model works + no formula labels", {
  skip_on_cran()
    # load data + fit model
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  fit <- mjoint(
    formLongFixed = list(log.b ~ year),
    formLongRandom = list(~ 1 | id),
    formSurv = Surv(years, status2) ~ age + edema,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 0.5, burnin = 20),
    verbose = FALSE)
  # tests
  expect_is(fit, "mjoint")
  expect_true(fit$conv)
  expect_equal(length(fixef(fit)), 2)
  expect_equal(nrow(ranef(fit)), fit$dims$n)
  expect_output(str(coef(fit)), "List of 5")
  expect_output(print(fit))
  expect_output(print(summary(fit)))
  expect_equal(length(fixef(fit, process = "Event")), 4)
})

test_that("multivariate model", {
  skip_on_cran()
  skip_on_os("mac")
    # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  set.seed(1)
  fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    control = list(convCrit = "sas", rav = 0.05, burnin = 20),
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


test_that("different convergence criteria", {
  skip_on_cran()
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
    control = list(convCrit = "sas", burnin = 4, mcmaxIter = 6),
    verbose = FALSE)
  fit2 <- update(fit1, control = list(convCrit = "rel", burnin = 4,
                                      mcmaxIter = 6))
  fit3 <- update(fit1, control = list(convCrit = "abs", burnin = 4,
                                      mcmaxIter = 6))
  fit4 <- update(fit1, control = list(convCrit = "either", burnin = 4,
                                      mcmaxIter = 6))
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
  skip_on_cran()
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
    control = list(convCrit = "abs", burnin = 4, mcmaxIter = 6, tol0 = 0.1),
    verbose = FALSE)
  expect_false(nrow(hvd1) == nrow(hvd2))
  expect_is(fit, "mjoint")
})


test_that("Gauss-Newton updates", {
  skip_on_cran()
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  set.seed(1)
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
                   burnin = 5, mcmaxIter = 7, gammaOpt = "GN"),
    verbose = FALSE)
  # tests
  expect_is(fit, "mjoint")
  expect_is(update(fit,
                   formSurv = Surv(fuyrs, status) ~ 1,
                   inits = list("gamma" = c(1.4, 0.8))),
            "mjoint")
})


test_that("no covariates in survival model", {
  skip_on_cran()
    # load data + fit model
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  fit <- mjoint(
    formLongFixed = list(log.b ~ year),
    formLongRandom = list(~ 1 | id),
    formSurv = Surv(years, status2) ~ 1,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 0.1, burnin = 3),
    verbose = FALSE)
  # tests
  expect_is(fit, "mjoint")
  expect_warning(baseHaz(fit), "No covariates in model to centre.")
})
