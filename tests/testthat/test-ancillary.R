library(joineRML)
context("Ancillary functions")

test_that("simulated data: intslope", {
  # simulate data
  beta <- rbind(c(0.5, 2, 1, 1),
                c(2, 2, -0.5, -1))
  D <- diag(4)
  D[1, 1] <- D[3, 3] <- 0.5
  D[1, 2] <- D[2, 1] <- D[3, 4] <- D[4, 3] <- 0.1
  D[1, 3] <- D[3, 1] <- 0.01
  sim <- simData(n = 250, beta = beta, D = D, sigma2 = c(0.25, 0.25),
                 censlam = exp(-0.2), gamma.y = c(-0.2, 1), ntms = 8)
  # tests
  expect_is(sim, "list")
  expect_output(str(sim), "List of 2")
  expect_equal(names(sim), c("longdat", "survdat"))
  expect_equal(ncol(sim$longdat), 8)
  expect_equal(dim(sim$survdat), c(250, 5))
})


test_that("simulated data: int", {
  # simulate data
  beta <- rbind(c(0.5, 2, 1, 1),
                c(2, 2, -0.5, -1))
  sim <- simData(n = 250, beta = beta, sigma2 = c(0.25, 0.25),
                 censlam = exp(-0.2), gamma.y = c(-0.2, 1), ntms = 8,
                 model = "int")
  # tests
  expect_is(sim, "list")
  expect_output(str(sim), "List of 2")
  expect_equal(names(sim), c("longdat", "survdat"))
  expect_equal(ncol(sim$longdat), 8)
  expect_equal(dim(sim$survdat), c(250, 5))
})


test_that("convergence plots", {
  # load data + fit model
  set.seed(1)
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  fit <- mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 5e-01, burnin = 5),
    verbose = FALSE)
  # tests
  expect_silent(plotConvergence(fit, params = "gamma"))
  expect_silent(plot(fit, type = "convergence", params = "beta"))
  expect_silent(plot(fit, type = "convergence", params = "gamma"))
  expect_silent(plot(fit, type = "convergence", params = "D"))
  expect_silent(plot(fit, type = "convergence", params = "sigma2"))
  expect_silent(plot(fit, type = "convergence", params = "loglik", discard = TRUE))
})


test_that("ranef plots + sampling", {
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  set.seed(1)
  fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
                 formLongRandom = ~ time | num,
                 formSurv = Surv(fuyrs, status) ~ age,
                 data = hvd,
                 timeVar = "time",
                 control = list(burnin = 6, tol0 = 5e-01))
  p <- plot(ranef(fit1, postVar = TRUE))
  # tests
  expect_true(is.ggplot(p))
  expect_error(sampleData(fit1, size = 1000, replace = FALSE), "Cannot select more subjects than in data without replacement")
})


test_that("dynamic predictions, residuals, fitted values, baseline hazard", {
  # Takes the most time so skip testing to pass Windows time-limit threshold
  skip_on_cran()
  skip_on_os("mac")
  # load data + fit model
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  set.seed(1)
  fit2 <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    control = list("burnin" = 10, mcmaxIter = 120, tol0 = 0.5),
    timeVar = "time")
  hvd2 <- droplevels(hvd[hvd$num == 1, ])
  test1 <- dynLong(fit2, hvd2)
  test2 <- dynLong(fit2, hvd2, u = 7)
  test3 <- dynSurv(fit2, hvd2)
  test4 <- dynSurv(fit2, hvd2, u = 7)
  test5 <- dynLong(fit2, hvd2, type = "simulated", M = 3)
  test6 <- dynSurv(fit2, hvd2, type = "simulated", M = 3)
  # tests: dynamic predictions
  expect_is(test1, "dynLong")
  expect_output(str(test1$pred), "List of 2")
  expect_silent(plot(test1))
  expect_output(print(test1))
  expect_is(test2, "dynLong")
  expect_is(test3, "dynSurv")
  expect_output(str(test3$pred), "data.frame")
  expect_silent(plot(test3))
  expect_output(print(test3))
  expect_is(test4, "dynSurv")
  expect_is(test5, "dynLong")
  expect_is(test6, "dynSurv")
  # tests: residuals + fitted values
  expect_output(str(resid(fit2, level = 0)), "List of 2")
  expect_output(str(resid(fit2, level = 1)), "List of 2")
  expect_output(str(fitted(fit2, level = 0)), "List of 2")
  expect_output(str(fitted(fit2, level = 1)), "List of 2")
  expect_error(fitted(fit2, level = 3))
  expect_equal(names(resid(fit2)), c("grad", "lvmi"))
  # tests: baseline hazard
  expect_is(baseHaz(fit2, se = TRUE), "data.frame")
  expect_is(baseHaz(fit2, centered = FALSE), "data.frame")
  # tests: missingg damts
  fit2$dmats <- NULL
  expect_error(fitted(fit2))
})

