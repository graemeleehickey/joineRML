test_that("SAS convergence treats increases and decreases symmetrically", {
  old <- list(D = matrix(1), beta = 1, sigma2 = 1, gamma = 0, haz = 0.1)
  con <- list(tol0 = 0.01, tol1 = 0.001, tol2 = 0.05, rav = 0.1, convCrit = "sas")
  check <- function(previous, next.values, verbose = FALSE) {
    joineRML:::convMonitor(previous, next.values, -100, -100, con, verbose)$conv
  }
  for (change in c(-100, -0.011, -0.01, 0.01, 0.011, 100)) {
    new <- old; new$gamma <- change
    expect_false(check(old, new))
  }
  for (change in c(-0.009, 0, 0.009)) {
    new <- old; new$gamma <- change
    expect_true(check(old, new))
  }
  old$gamma <- con$rav
  new <- old; new$gamma <- old$gamma + 0.007
  expect_true(check(old, new))
  old$gamma <- con$rav + 0.001
  new <- old; new$gamma <- old$gamma + 0.007
  expect_false(check(old, new))
  old$gamma <- 1
  for (direction in c(-1, 1)) {
    new <- old; new$gamma <- 1 + direction * 0.04
    expect_true(check(old, new))
    new$gamma <- 1 + direction * 0.06
    expect_false(check(old, new))
  }

  old$gamma <- 0
  new <- old; new$gamma <- 0.02
  reordered <- old[c("haz", "gamma", "sigma2", "D", "beta")]
  expect_identical(check(reordered, new), check(old, new))
  output <- capture.output(check(old, new, verbose = TRUE))
  small.change <- grep("for parameters  <", output, value = TRUE)
  expect_match(small.change, "= 0.02")
})

test_that("convergence still requires burn-in and three successive checks", {
  for (criterion in c("sas", "rel")) {
    fit <- suppressMessages(mjoint(
      y1 ~ time, ~ 1 | id, Surv(T, status) ~ 1,
      data = small_joint_data(), timeVar = "time", inits = list(gamma = 0.1),
      control = list(type = "halton", nMC = 128, nMCmax = 128, burnin = 2,
                     mcmaxIter = 8, convCrit = criterion, tol0 = Inf, tol2 = Inf)
    ))
    expect_true(fit$conv)
    expect_equal(ncol(fit$history), 5L)
  }
})
