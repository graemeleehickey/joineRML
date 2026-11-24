test_that("counts describe subjects, events, and each measurement schedule", {
  dat <- small_joint_data()
  fit <- fit_small_joint(
    data = list(dat, subset(dat, !(id %% 2 == 0 & time == 1))),
    formLongFixed = list(first = y1 ~ time, second = y2 ~ time),
    formLongRandom = list(first = ~ 1 | id, second = ~ 1 | id),
    pfs = FALSE, mcmaxIter = 2)
  expect_equal(nobs(fit), 30)
  expect_equal(nobs(fit), attr(logLik(fit), "nobs"))
  g <- glance(fit)
  expect_equal(g$nobs, 30)
  expect_equal(g$nevent, 20)
  expect_equal(g$nobs_longitudinal_1, 90)
  expect_equal(g$nobs_longitudinal_2, 75)
  expect_equal(g$npar, 11)
  expect_false(g$converged)
  expect_equal(g$AIC, AIC(fit))
  expect_equal(g$BIC, BIC(fit))
  expect_equal(g$logLik, as.numeric(logLik(fit)))
  expect_error(augment(fit), "identical|different")
})

test_that("coefficient and covariance selections agree and preserve defaults", {
  fit <- fit_small_joint(
    formLongFixed = list(first = y1 ~ time, second = y2 ~ time),
    formLongRandom = list(first = ~ 1 | id, second = ~ 1 | id))
  expect_identical(coef(fit), fit$coefficients)
  expect_identical(coef(fit, complete = FALSE), fit$coefficients)
  expect_identical(coef(fit, process = "Longitudinal"), fit$coefficients$beta)
  expect_identical(coef(fit, process = "Event"), fit$coefficients$gamma)
  full <- vcov(fit)
  for (k in 1:2) {
    cf <- coef(fit, process = "Longitudinal", k = k)
    v <- vcov(fit, process = "Longitudinal", k = k)
    expect_named(cf, colnames(model.matrix(fit, k = k)))
    expect_identical(rownames(v), names(cf))
    expect_identical(colnames(v), names(cf))
    expect_equal(unname(v), unname(full[3 + (2 * k - 1:0), 3 + (2 * k - 1:0)]))
    expect_equal(as.vector(model.matrix(fit, k = k) %*% cf),
                 as.vector(fitted(fit)[[k]]))
    expect_equal(vcov(fit, correlation = TRUE, process = "Longitudinal", k = k),
                 cov2cor(v))
  }
  ev <- vcov(fit, process = "Event")
  expect_identical(rownames(ev), names(coef(fit, process = "Event")))
  expect_equal(unname(ev), unname(full[10:11, 10:11]))
  expect_equal(dim(model.matrix(fit, process = "Event")), c(30L, 0L))
  single <- fit_small_joint(formLongFixed = list(y1 = y1 ~ 1))
  expect_equal(dim(vcov(single, process = "Longitudinal", k = 1)), c(1L, 1L))
  expect_equal(dim(vcov(single, process = "Event")), c(1L, 1L))
  expect_error(coef(fit, k = 1), "Specify")
  expect_error(vcov(fit, k = 1), "Specify")
  expect_error(coef(fit, process = "Event", k = 1), "only used")
  for (k in list(0, 3, 1.5, NA, Inf, "first", c(1, 2))) {
    expect_error(model.frame(fit, k = k), "'k'")
    expect_error(model.matrix(fit, k = k), "'k'")
    expect_error(terms(fit, k = k), "'k'")
    expect_error(formula(fit, k = k), "'k'")
    expect_error(coef(fit, process = "Longitudinal", k = k), "'k'")
  }
})

test_that("frames and matrices retain transformations and fitted contrasts", {
  dat <- small_joint_data()
  dat$group <- factor(rep(c("a", "b", "c"), each = 2, length.out = 30)[dat$id])
  dat$age <- dat$id + 30
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  fit <- fit_small_joint(data = dat, pfs = FALSE,
    formLongFixed = list(y1 = y1 ~ poly(time, 2) * group),
    formSurv = Surv(T, status) ~ group + I(age / 10))
  frame <- model.frame(fit)
  matrix <- model.matrix(fit)
  event <- model.matrix(fit, process = "Event")
  expect_equal(unname(model.response(frame)), fit$data[[1]]$y1)
  expect_equal(model.response(model.frame(fit, process = "Event")), fit$sfit$y,
               ignore_attr = TRUE)
  expect_equal(unname(matrix), unname(do.call(rbind, fit$dmats$l$Xik.list[[1]])),
               ignore_attr = TRUE)
  expect_equal(event, fit$sfit$x)
  expect_true(all(colnames(event) %in% names(coef(fit, process = "Event"))))
  expect_s3_class(terms(fit), "terms")
  expect_equal(attr(terms(fit), "term.labels"),
               c("poly(time, 2)", "group", "poly(time, 2):group"))
  expect_s3_class(terms(fit, process = "Event"), "terms")
  options(contrasts = c("contr.treatment", "contr.poly"))
  dat$time <- dat$time + 100
  expect_identical(model.frame(fit), frame)
  expect_identical(model.matrix(fit), matrix)
  expect_identical(model.matrix(fit, process = "Event"), event)
  older <- fit
  older$model.frames <- NULL
  expect_equal(model.matrix(older), matrix)
  expect_equal(model.frame(older, process = "Event"),
               model.frame(fit, process = "Event"), ignore_attr = TRUE)
  expect_error(model.matrix(fit, data = dat), "Unused argument")
  expect_error(vcov(fit), "No approximate covariance")
})
