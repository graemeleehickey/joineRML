library(joineRML)

test_that("baseHaz returns centered and uncentered estimates", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()

  bh <- baseHaz(fit)
  expect_named(bh, c("time", "haz"))
  expect_equal(nrow(bh), length(fit$dmats$t$tj))

  bh_unc <- baseHaz(fit, centered = FALSE)
  expect_named(bh_unc, c("time", "haz"))
  expect_false(isTRUE(all.equal(bh$haz, bh_unc$haz)))

  bh_se <- baseHaz(fit, se = TRUE)
  expect_named(bh_se, c("time", "haz", "se"))
  expect_true(all(bh_se$se >= 0))
})

test_that("baseHaz input checks", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  expect_error(baseHaz(1), "Use only with 'mjoint' model objects")
  expect_error(baseHaz(fit, centered = FALSE, se = TRUE),
               "Can only estimate standard errors for the centered case")
})

test_that("confint returns expected dimensions per parameter set", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()

  expect_equal(ncol(confint(fit)), 2)
  expect_equal(nrow(confint(fit, parm = "Longitudinal")), fit$dims$p)
  expect_equal(nrow(confint(fit, parm = "Event")),
               with(fit$dims, q + K))
})

test_that("accessor and print methods run", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()

  expect_type(sigma(fit), "double")
  expect_length(fixef(fit, process = "Longitudinal"), fit$dims$p)
  expect_length(fixef(fit, process = "Event"), with(fit$dims, q + K))
  expect_s3_class(formula(fit, process = "Longitudinal"), "formula")
  expect_output(print(fit))
  expect_output(print(summary(fit)))
})

test_that("print.summary.mjoint covers univariate branch and input check", {
  skip_on_cran()
  skip_on_os("mac")
  fit <- fit_pbc_uni()
  fit.summ <- summary(fit)

  expect_s3_class(fit.summ, "summary.mjoint")
  expect_output(print(fit.summ), "Univariate linear mixed-effects model")
  expect_output(print(fit.summ), "Cox proportional hazards model")
  expect_output(print(fit.summ), "method: approx")
  expect_error(print.summary.mjoint(1), "Use only with 'summary.mjoint' objects")
})

