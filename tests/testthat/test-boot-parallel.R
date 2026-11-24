test_that("parallel bootstrap preserves existing settings and closes its connections", {
  skip_on_cran()
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")
  # The user's existing foreach setup must remain usable after either outcome.
  doParallel::registerDoParallel(cores = 1L)
  on.exit({
    doParallel::stopImplicitCluster()
    foreach::registerDoSEQ()
  }, add = TRUE)
  settings <- function() list(name = foreach::getDoParName(),
                              workers = foreach::getDoParWorkers())
  sockets <- function() sum(showConnections(all = TRUE)[, "class"] == "sockconn")
  before <- settings()
  connections <- sockets()
  fit <- fit_small_joint()

  bs <- bootSE(fit, nboot = 2, ncores = 2, progress = FALSE)
  expect_equal(bs$conv, 2)
  expect_true(all(is.finite(bs$beta.se)))
  expect_equal(settings(), before)
  expect_equal(sockets(), connections)
  expect_identical(foreach::`%dopar%`(foreach::foreach(i = 1:2), i + 1L),
                    list(2L, 3L))

  # Force a fitting error in a worker, after the new connections are opened.
  fit$formSurv <- Surv(T, status) ~ missing_covariate
  expect_error(bootSE(fit, nboot = 2, ncores = 2, progress = FALSE),
                "missing_covariate")
  expect_equal(settings(), before)
  expect_equal(sockets(), connections)
  expect_identical(foreach::`%dopar%`(foreach::foreach(i = 1:2), i + 1L),
                    list(2L, 3L))
})
