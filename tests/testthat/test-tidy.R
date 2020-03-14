# test tidy, augment, glance methods for mjoint object (joineRML package)
context("Tidiers")
library(joineRML)

# Data
data(heart.valve, package = "joineRML")
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi) & heart.valve$num <= 50, ]

# Model fits
fit1 <- suppressMessages(joineRML::mjoint(
  formLongFixed = list("grad" = log.grad ~ time),
  formLongRandom = list("grad" = ~ 1 | num),
  formSurv = survival::Surv(fuyrs, status) ~ age,
  data = hvd,
  timeVar = "time"
))
fit2 <- suppressMessages(joineRML::mjoint(
  formLongFixed = list(
    "grad" = log.grad ~ time,
    "lvmi" = log.lvmi ~ time
  ),
  formLongRandom = list(
    "grad" = ~ 1 | num,
    "lvmi" = ~ 1 | num
  ),
  formSurv = Surv(fuyrs, status) ~ age,
  data = hvd,
  timeVar = "time"
))

# Bootstrapped SEs
set.seed(092345798)
bSE1 <- suppressMessages(bootSE(fit1, nboot = 3, safe.boot = TRUE, progress = FALSE))
bSE2 <- suppressMessages(bootSE(fit2, nboot = 3, safe.boot = TRUE, progress = FALSE))

test_that("tidy works on mjoint models with a single longitudinal process", {
  td <- tidy(fit1)
  expect_equal(nrow(td), 2)
  td <- tidy(fit1, component = "survival")
  expect_equal(nrow(td), 2)
  td <- tidy(fit1, component = "longitudinal")
  expect_equal(nrow(td), 2)
  td <- tidy(fit1, component = "survival", bootSE = bSE1)
  expect_equal(nrow(td), 2)
  td <- tidy(fit1, component = "longitudinal", bootSE = bSE1)
  expect_equal(nrow(td), 2)

  td <- tidy(fit1, component = "survival")
  expect_equal(td$term, c("age", "gamma_1"))
  td <- tidy(fit1, component = "longitudinal")
  expect_equal(td$term, c("(Intercept)_1", "time_1"))
})

test_that("tidy works on mjoint models with more than one longitudinal process", {
  td <- tidy(fit2)
  expect_equal(nrow(td), 3)
  td <- tidy(fit2, component = "survival")
  expect_equal(nrow(td), 3)
  td <- tidy(fit2, component = "longitudinal")
  expect_equal(nrow(td), 4)
  td <- tidy(fit2, component = "survival", bootSE = bSE2)
  expect_equal(nrow(td), 3)
  td <- tidy(fit2, component = "longitudinal", bootSE = bSE2)
  expect_equal(nrow(td), 4)

  td <- tidy(fit2, component = "survival")
  expect_equal(td$term, c("age", "gamma_1", "gamma_2"))
  td <- tidy(fit2, component = "longitudinal")
  expect_equal(td$term, c("(Intercept)_1", "time_1", "(Intercept)_2", "time_2"))
})

test_that("augment works on mjoint models with a single longitudinal process", {
  augdf <- augment(fit1)
  expect_equal(nrow(augdf), nrow(hvd))
  expect_equal(ncol(augdf), ncol(hvd) + 4)
  expect_equal(names(augdf), c(names(hvd), ".fitted_grad_0", ".fitted_grad_1", ".resid_grad_0", ".resid_grad_1"))
})

test_that("augment works on mjoint models with more than one longitudinal process", {
  augdf <- augment(fit2)
  expect_equal(nrow(augdf), nrow(hvd))
  expect_equal(ncol(augdf), ncol(hvd) + 8)
  expect_equal(names(augdf), c(names(hvd), ".fitted_grad_0", ".fitted_lvmi_0", ".fitted_grad_1", ".fitted_lvmi_1", ".resid_grad_0", ".resid_lvmi_0", ".resid_grad_1", ".resid_lvmi_1"))
})

test_that("augment returns the same output whether we pass 'data' or not", {
  expect_equal(object = names(augment(fit1)), expected = names(augment(fit1, data = list(hvd))))
  expect_equal(object = dim(augment(fit1)), expected = dim(augment(fit1, data = list(hvd))))
  expect_equal(object = names(augment(fit2)), expected = names(augment(fit2, data = list(hvd))))
  expect_equal(object = dim(augment(fit2)), expected = dim(augment(fit2, data = list(hvd))))
})

test_that("glance works on mjoint models with a single longitudinal process", {
  glnc <- glance(fit1)
  expect_equal(nrow(glnc), 1)
  expect_equal(ncol(glnc), 4)
  expect_equal(names(glnc), c("sigma2_1", "AIC", "BIC", "logLik"))
})

test_that("glance works on mjoint models with more than one longitudinal process", {
  glnc <- glance(fit2)
  expect_equal(nrow(glnc), 1)
  expect_equal(ncol(glnc), 5)
  expect_equal(names(glnc), c("sigma2_1", "sigma2_2", "AIC", "BIC", "logLik"))
})

test_that("tidy fails if passing a bootSE object that is not a bootSE object", {
  expect_error(tidy(fit1, bootSE = list()))
  expect_error(tidy(fit2, bootSE = list()))
})

test_that("tidy returns confidence intervals if required", {
  expect_true(all(c("conf.low", "conf.high") %in% names(tidy(fit1, conf.int = TRUE))))
  expect_true(all(c("conf.low", "conf.high") %in% names(tidy(fit1, conf.int = TRUE, bootSE = bSE1))))
  expect_true(all(c("conf.low", "conf.high") %in% names(tidy(fit2, conf.int = TRUE))))
  expect_true(all(c("conf.low", "conf.high") %in% names(tidy(fit2, conf.int = TRUE, bootSE = bSE2))))
})

test_that("augment fails if cannot extract data from x", {
  fit1_broken <- fit1
  fit1_broken$data <- NULL
  expect_error(augment(fit1_broken))
  nd <- list(fit1$data[[1]], fit1$data[[1]][sample(x = seq(nrow(fit1$data[[1]])), size = nrow(fit1$data[[1]])), ])
  expect_error(augment(fit1_broken))
  fit2_broken <- fit2
  fit2_broken$data <- NULL
  expect_error(augment(fit2_broken))
  nd <- list(fit2$data[[1]], fit2$data[[1]][sample(x = seq(nrow(fit2$data[[1]])), size = nrow(fit2$data[[1]])), ])
  expect_error(augment(fit2_broken))
})
