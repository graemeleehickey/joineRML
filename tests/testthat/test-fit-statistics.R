# Independently integrate the random-intercept model with the same fixed nodes.
random_intercept_loglik <- function(fit) {
  theta <- coef(fit)
  dat <- fit$data[[1]]
  surv <- fit$survData
  times <- sort(unique(surv$T[surv$status == 1]))
  nodes <- randtoolbox::halton(fit$finalnMC, dim = 1, normal = TRUE)
  sum(vapply(seq_len(nrow(surv)), function(i) {
    d <- dat[dat$id == surv$id[i], ]
    mu <- as.vector(model.matrix(fit$formLongFixed[[1]], d) %*% theta$beta)
    m <- nrow(d)
    variance <- theta$D[1, 1]
    sigma2 <- theta$sigma2[1]
    log.long <- mvtnorm::dmvnorm(d$y1, mu,
      sigma = diag(sigma2, m) + matrix(variance, m, m), log = TRUE)
    post.var <- 1 / (1 / variance + m / sigma2)
    post.mean <- post.var * sum(d$y1 - mu) / sigma2
    b <- post.mean + sqrt(post.var) * as.vector(nodes)
    eta <- theta$gamma[1] * b
    H <- sum(theta$haz[times <= surv$T[i]])
    log.surv <- -H * exp(eta)
    if (surv$status[i] == 1) {
      log.surv <- log.surv + log(theta$haz[match(surv$T[i], times)]) + eta
    }
    log.long + log(mean(exp(log.surv)))
  }, numeric(1)))
}

test_that("reported likelihood uses the returned coefficients", {
  for (pfs in c(FALSE, TRUE)) {
    for (iterations in c(2, 4)) {
      fit <- fit_small_joint(pfs = pfs, mcmaxIter = iterations)
      expect_identical(fit$conv, iterations == 4)
      expect_equal(as.numeric(logLik(fit)), random_intercept_loglik(fit),
                   tolerance = 1e-10)
      expect_gt(abs(fit$log.lik - fit$ll.hx[iterations]), 1e-6)
      expect_equal(unname(tail(fit$history["long_time_1", ], 1)),
                   unname(coef(fit)$beta["time_1"]))
    }
  }
})

test_that("likelihood-only evaluation does not require a parameter update", {
  fit <- fit_small_joint()
  l <- fit$dmats$l
  l$XtX.inv <- NULL
  out <- joineRML:::stepEM(coef(fit), l, fit$dmats$t, fit$dmats$z,
    nMC = fit$finalnMC, verbose = FALSE, gammaOpt = "NR", pfs = FALSE,
    type = "halton", likelihood.only = TRUE)
  expect_named(out, "ll")
  expect_equal(out$ll, random_intercept_loglik(fit), tolerance = 1e-10)
})

test_that("BIC consistently uses the number of subjects", {
  fit <- fit_small_joint()
  ll <- logLik(fit)
  expected <- -2 * as.numeric(ll) + log(fit$dims$n) * attr(ll, "df")
  expect_equal(attr(ll, "nobs"), fit$dims$n)
  expect_equal(BIC(fit), expected)
  expect_equal(summary(fit)$BIC, expected)
  expect_equal(glance(fit)$BIC, expected)
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
})
