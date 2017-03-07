library(joineRML)
context("Ancillary functions")

test_that("simulated data", {
  # simulate data
  beta <- rbind(c(0.5, 2, 1, 1),
                c(2, 2, -0.5, -1))
  D <- diag(4)
  D[1, 1] <- D[3, 3] <- 0.5
  D[1, 2] <- D[2, 1] <- D[3, 4] <- D[4, 3] <- 0.1
  D[1, 3] <- D[3, 1] <- 0.01
  sim <- simData(n = 250, beta = beta, D = D, sigma2 = c(0.25, 0.25),
                 censlam = exp(-0.2), gamma.y = c(-.2, 1), ntms = 8)
  # tests
  expect_is(sim, "list")
  expect_output(str(sim), "List of 2")
  expect_equal(names(sim), c("longdat", "survdat"))
  expect_equal(ncol(sim$longdat), 8)
  expect_equal(dim(sim$survdat), c(250, 5))
})


test_that("convergence plots", {
  # load data + fit model
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  fit <- mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 1e-02),
    verbose = FALSE)
  # tests
  expect_silent(plotConvergence(fit, params = "gamma"))
  expect_silent(plot(fit, type = "convergence", params = "beta"))
  expect_silent(plot(fit, type = "convergence", params = "gamma"))
  expect_silent(plot(fit, type = "convergence", params = "D"))
  expect_silent(plot(fit, type = "convergence", params = "sigma2"))
})
