test_that("unsupported survival terms are rejected before fitting", {
  dat <- small_joint_data()
  dat$x <- dat$id / 30
  dat$group <- factor(dat$id %% 2)
  fit <- function(formula, separate = FALSE) {
    mjoint(y1 ~ time, ~ 1 | id, formula, data = dat,
           survData = if (separate) dat[!duplicated(dat$id), ] else NULL,
           timeVar = "time")
  }
  terms <- c("strata(group)", "offset(x)", "cluster(id)", "tt(x)",
             "frailty(group)", "pspline(x)", "ridge(x)",
             "frailty.gamma(group)", "frailty.gaussian(group)", "frailty.t(group)")
  for (term in terms) {
    for (separate in c(FALSE, TRUE)) {
      formula <- as.formula(paste("Surv(T, status) ~", term))
      expect_error(fit(formula, separate), "Unsupported survival formula term")
    }
  }
  for (term in c("survival::strata(group)", "stats::offset(x)",
                 "survival::cluster(id)", "survival::tt(x)",
                 "survival::frailty(group)", "survival::pspline(x)",
                 "survival::ridge(x)", "survival:::strata(group)",
                 "I(offset(x))", "x:strata(group)", "x * strata(group)")) {
    expect_error(fit(as.formula(paste("survival::Surv(T, status) ~", term))),
                  "Unsupported survival formula term")
  }
  expect_error(fit(Surv(T, status) ~ x + offset(x) + strata(group)),
                "offset\\(\\), strata\\(\\)")
  # tt() does not need to be defined, nor its argument evaluated, for rejection.
  expect_error(fit(Surv(T, status) ~ tt(missing_covariate)),
                "Unsupported survival formula term.*tt")
})

test_that("penalised survival terms cannot bypass the check through an alias", {
  dat <- small_joint_data()
  dat$x <- dat$id / 30
  custom_penalty <- function(x) survival::ridge(x)
  expect_error(mjoint(y1 ~ time, ~ 1 | id,
                       Surv(T, status) ~ custom_penalty(x), dat, timeVar = "time"),
                "Penalised survival terms.*custom_penalty")
})

test_that("ordinary baseline covariates and transformations remain supported", {
  dat <- small_joint_data()
  dat$x <- dat$id / 30
  dat$group <- factor(dat$id %% 2)
  dat$strata <- sin(dat$id)
  dat$offset <- cos(dat$id)
  dat$X <- cbind(dat$x, dat$x^2)
  for (formula in list(Surv(T, status) ~ 1,
                       survival::Surv(T, status) ~ x * group + I(x^2),
                       Surv(T, status) ~ splines::ns(x, df = 2),
                       Surv(T, status) ~ strata + offset,
                       Surv(T, status) ~ X[, 1])) {
    if (any(attr(terms(formula), "order") > 1L)) {
      # The preliminary Cox survival curve warns for interaction terms.
      expect_warning(fit <- fit_small_joint(dat, formSurv = formula),
                      "model contains interactions")
    } else {
      fit <- fit_small_joint(dat, formSurv = formula)
    }
    expect_true(fit$conv)
    expect_true(all(is.finite(unlist(coef(fit)))))
    expect_true(is.finite(as.numeric(logLik(fit))))
    expected <- model.matrix(formula, fit$survData)[, -1L, drop = FALSE]
    expect_equal(unname(model.matrix(fit, process = "Event")), unname(expected),
                  ignore_attr = TRUE)
  }
})
