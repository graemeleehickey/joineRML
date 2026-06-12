library(joineRML)

quick_fit <- function() {
  data(pbc2)
  pbc2$log.b <- log(pbc2$serBilir)
  set.seed(4821)
  mjoint(
    formLongFixed = list("log.b" = log.b ~ year),
    formLongRandom = list("log.b" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year",
    control = list(convCrit = "abs", tol0 = 0.5, burnin = 20, mcmaxIter = 100),
    verbose = FALSE)
}

test_that("plot.dynSurv runs for first-order and simulated predictions", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- quick_fit()
  nd <- droplevels(subset(pbc2, id == 1))
  ds <- dynSurv(fit, nd, progress = FALSE)
  set.seed(11)
  ds_sim <- dynSurv(fit, nd, type = "simulated", M = 5, progress = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(suppressWarnings(plot(ds, main = "P1")))
  expect_invisible(suppressWarnings(plot(ds, smooth = TRUE)))
  expect_invisible(suppressWarnings(plot(ds_sim)))
  expect_invisible(suppressWarnings(plot(ds_sim, estimator = "mean")))
  expect_error(suppressWarnings(plot(ds_sim, estimator = "bad")),
               "'mean' or 'median'")
  expect_error(suppressWarnings(plot(ds, ylab1 = c("a", "b"))),
               "does not match number of outcomes")
  expect_error(plot.dynSurv(1), "Use only with 'dynSurv' objects")
})

test_that("plot.dynLong runs for first-order and simulated predictions", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- quick_fit()
  nd <- droplevels(subset(pbc2, id == 1))
  dl <- dynLong(fit, nd, progress = FALSE)
  set.seed(11)
  dl_sim <- dynLong(fit, nd, type = "simulated", M = 5, progress = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(dl, main = "P1"))
  expect_invisible(plot(dl_sim))
  expect_invisible(plot(dl_sim, estimator = "mean"))
  expect_error(plot(dl, ylab = c("a", "b")),
               "does not match number of outcomes")
  expect_error(plot.dynLong(1), "Use only with 'dynLong' objects")
})
