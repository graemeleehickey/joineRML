library(joineRML)

test_that("dynSurv first-order prediction", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  nd <- pbc_newdata()

  ds <- dynSurv(fit, nd, progress = FALSE)
  expect_s3_class(ds, "dynSurv")
  expect_equal(ds$type, "first-order")
  expect_named(ds$pred, c("u", "surv"))
  expect_true(all(ds$pred$surv >= 0 & ds$pred$surv <= 1))
  expect_output(print(ds))

  ds_u <- dynSurv(fit, nd, u = 10, progress = FALSE)
  expect_equal(ds_u$pred$u, 10)

  ds_h <- dynSurv(fit, nd, horizon = 2, progress = FALSE)
  expect_equal(ds_h$horizon, 2)
  expect_output(print(ds_h))
})

test_that("dynSurv simulated prediction", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  nd <- pbc_newdata()

  set.seed(11)
  ds <- dynSurv(fit, nd, type = "simulated", M = 5, progress = FALSE)
  expect_s3_class(ds, "dynSurv")
  expect_equal(ds$type, "simulated")
  expect_named(ds$pred, c("u", "mean", "median", "lower", "upper"))
  expect_equal(ds$M, 5)
  expect_true(ds$accept >= 0 && ds$accept <= 1)
  expect_output(print(ds))
})

test_that("dynSurv input checks", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  nd <- pbc_newdata()

  expect_error(dynSurv(fit, nd, u = 4, horizon = 2, progress = FALSE),
               "Cannot specify 'u' and 'horizon'")
  expect_error(dynSurv(fit, nd, u = -1, progress = FALSE),
               "must be greater than last observation time")
  expect_error(dynSurv(fit, nd, type = "nonsense", progress = FALSE),
               "first-order")
  expect_error(dynSurv(1, nd, progress = FALSE),
               "Use only with 'mjoint' model objects")
})

test_that("dynLong first-order and simulated prediction", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  nd <- pbc_newdata()

  dl <- dynLong(fit, nd, progress = FALSE)
  expect_s3_class(dl, "dynLong")
  expect_equal(dl$type, "first-order")
  expect_output(print(dl))

  dl_u <- dynLong(fit, nd, u = 10, progress = FALSE)
  expect_s3_class(dl_u, "dynLong")

  set.seed(11)
  dl_sim <- dynLong(fit, nd, type = "simulated", M = 5, progress = FALSE)
  expect_s3_class(dl_sim, "dynLong")
  expect_equal(dl_sim$type, "simulated")
  expect_output(print(dl_sim))
})
