# Use a known joint covariance to check the sampling approximation directly.
theta_draw_for_parameters <- function(theta, covariance, K = length(theta$sigma2)) {
  draw <- joineRML:::thetaDraw
  env <- new.env(parent = environment(draw))
  env$coef <- function(object) theta
  env$vcov <- function(object) covariance
  env$baseHaz <- function(object, se) data.frame(haz = 0.1, se = 0)
  environment(draw) <- env
  function() draw(list(dims = list(K = K, r = nrow(theta$D))))
}

theta_draw_parameters <- function() {
  list(D = matrix(0.8), beta = c(intercept = 2, time = 0.3),
       sigma2 = c(y1 = 0.4), gamma = c(y1 = 0.1), haz = 0.1)
}

test_that("joint draws have the intended log-variance moments and covariances", {
  theta <- theta_draw_parameters()
  covariance <- diag(c(0.01, 0.09, 0.04, 0.01, 0.02))
  covariance[4, c(2, 3, 5)] <- c(0.015, -0.006, 0.005)
  covariance[c(2, 3, 5), 4] <- covariance[4, c(2, 3, 5)]
  draw <- theta_draw_for_parameters(theta, covariance)
  set.seed(482)
  samples <- t(replicate(4000, {
    sample <- draw()
    c(sample$D, sample$beta, log(sample$sigma2), sample$gamma)
  }))

  # The derivatives and expected moments are specified independently.
  J <- diag(c(1, 1, 1, 2.5, 1))
  expected.mean <- c(0.8, 2, 0.3, log(0.4), 0.1)
  expected.covariance <- J %*% covariance %*% J
  n <- nrow(samples)
  mean.mcse <- sqrt(diag(expected.covariance) / n)
  covariance.mcse <- sqrt((expected.covariance^2 +
    outer(diag(expected.covariance), diag(expected.covariance))) / (n - 1))
  expect_true(all(abs(colMeans(samples) - expected.mean) < 6 * mean.mcse))
  expect_true(all(abs(cov(samples) - expected.covariance) < 6 * covariance.mcse))
})

test_that("multivariate draws preserve parameter order, names, and reproducibility", {
  theta <- list(D = matrix(c(1, 0.2, 0.2, 2), 2),
                beta = c(y1.intercept = 2, y1.time = 0.3, y2.intercept = 1),
                sigma2 = c(y1 = 0.4, y2 = 0.7),
                gamma = c(y1 = 0.1, y2 = -0.1), haz = 0.1)
  covariance <- diag(0.002, 10)
  draw <- theta_draw_for_parameters(theta, covariance)
  set.seed(231)
  first <- draw()
  set.seed(231)
  expect_identical(draw(), first)
  expect_named(first, names(theta))
  for (component in c("beta", "sigma2", "gamma")) {
    expect_named(first[[component]], names(theta[[component]]))
  }
  expect_equal(dim(first$D), c(2L, 2L))
  expect_true(all(is.finite(first$sigma2) & first$sigma2 > 0))

  fixed.draw <- theta_draw_for_parameters(theta, matrix(0, 10, 10))
  expect_equal(fixed.draw(), theta, tolerance = 1e-12)
})

test_that("invalid fitted variances and covariance matrices give clear errors", {
  theta <- theta_draw_parameters()
  for (variance in c(0, -0.1, NA, NaN, Inf)) {
    invalid <- theta
    invalid$sigma2[] <- variance
    draw <- theta_draw_for_parameters(invalid, diag(5))
    expect_error(draw(), "finite, strictly positive fitted residual variances")
  }
  draw <- theta_draw_for_parameters(theta, diag(5), K = 2)
  expect_error(draw(), "fitted residual variances")
  invalid <- theta
  invalid$beta[1] <- Inf
  expect_error(theta_draw_for_parameters(invalid, diag(5))(),
               "finite fitted parameters")

  asymmetric <- diag(5)
  asymmetric[1, 2] <- 0.5
  for (covariance in list(diag(4), rep(1, 25), matrix(NA_real_, 5, 5),
                          matrix(Inf, 5, 5), asymmetric)) {
    expect_error(theta_draw_for_parameters(theta, covariance)(),
                 "finite, symmetric covariance matrix matching")
  }
  indefinite <- diag(5)
  indefinite[1, 2] <- indefinite[2, 1] <- 2
  expect_error(theta_draw_for_parameters(theta, indefinite)(),
               "Cannot simulate parameters.*not positive semidefinite")

  invalid <- theta
  invalid$sigma2[] <- 1e-300
  expect_error(theta_draw_for_parameters(invalid, diag(5))(),
               "not finite after transforming residual variances")
})

test_that("overflow and underflow are reported instead of clipping variances", {
  covariance <- diag(5)
  covariance[4, 4] <- 1e300
  draw <- theta_draw_for_parameters(theta_draw_parameters(), covariance)
  for (seed in 1:2) {
    set.seed(seed)
    expect_error(draw(), "not finite or residual variances are not strictly positive")
  }
})

test_that("small-sample simulated predictions use positive residual variances", {
  fit <- fit_small_joint(subset(small_joint_data(), id <= 6))
  original.coefficients <- coef(fit)
  original.covariance <- vcov(fit)
  original.sigma <- sigma(fit)
  nd <- subset(fit$data[[1]], id == 1)

  # Only suppress the expected warning for an event model without covariates.
  without.centering.warning <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("No covariates in model to centre", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  set.seed(18)
  variances <- without.centering.warning(replicate(1000, joineRML:::thetaDraw(fit)$sigma2))
  expect_true(all(is.finite(variances) & variances > 0))
  set.seed(18)
  long <- without.centering.warning(dynLong(fit, nd, u = 3, type = "simulated",
                                          M = 1000, progress = FALSE))
  set.seed(18)
  surv <- without.centering.warning(dynSurv(fit, nd, u = 3.5, type = "simulated",
                                          M = 1000, progress = FALSE))
  expect_true(all(is.finite(as.matrix(long$pred[[1]]))))
  expect_true(all(is.finite(as.matrix(surv$pred))))
  probabilities <- as.matrix(surv$pred[, c("mean", "median", "lower", "upper")])
  expect_true(all(probabilities >= 0 & probabilities <= 1))
  expect_identical(coef(fit), original.coefficients)
  expect_identical(vcov(fit), original.covariance)
  expect_identical(sigma(fit), original.sigma)
})
