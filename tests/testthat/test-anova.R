# Small fixed-iteration fits test comparison calculations, not model adequacy.
comparison_models <- function() {
  dat <- small_joint_data()
  dat$group <- factor(rep(rep(c("A", "B"), 15), each = 3))
  dat$x <- rep(seq_len(30) / 30, each = 3)
  fit <- function(fixed, survival = Surv(T, status) ~ 1, pfs = TRUE) {
    suppressMessages(mjoint(fixed, ~ 1 | id, survival, data = dat,
      timeVar = "time", inits = list(gamma = rep(0.1, length(all.vars(survival)) - 1)),
      pfs = pfs, control = list(type = "halton", nMC = 128, nMCmax = 128,
        burnin = 1, mcmaxIter = 4, convCrit = "abs", tol0 = Inf)))
  }
  list(reduced = fit(y1 ~ 1), full = fit(y1 ~ time),
       extended = fit(y1 ~ time + group, Surv(T, status) ~ x, pfs = FALSE))
}

# Reuse the fits across checks in this file.
comparison_fixtures <- comparison_models()

test_that("nested comparisons report the expected likelihood-ratio calculations", {
  fits <- comparison_fixtures
  tab <- anova(fits$reduced, fits$full, nMC = 4096)
  expect_s3_class(tab, "anova")
  expect_s3_class(tab, "data.frame")
  expect_equal(tab$df, unname(vapply(fits[1:2], function(x) attr(logLik(x), "df"), numeric(1))))
  expect_equal(tab$Df, c(NA, 1))
  expect_gt(tab$L.Ratio[2], 0)
  expect_equal(tab$L.Ratio[2], 2 * diff(tab$logLik))
  expect_equal(tab[["Pr(>Chisq)"]][2], pchisq(tab$L.Ratio[2], 1, lower.tail = FALSE))
  expect_equal(tab$AIC, -2 * tab$logLik + 2 * tab$df)
  expect_equal(tab$BIC, -2 * tab$logLik + log(30) * tab$df)
  expect_equal(tab$LR.change[2], abs(tab$L.Ratio[2] -
    2 * diff(attr(tab, "integration")$logLik.small)))
  expect_equal(attr(tab, "integration")$nMC.small, 2048)
  expect_output(print(tab), "Halton")
})

test_that("descriptive comparisons use stored likelihoods and leave models unchanged", {
  fits <- comparison_fixtures
  original <- serialize(fits, NULL)
  tab <- anova(fits$full, fits$reduced, test = FALSE)
  expect_equal(tab$logLik, unname(vapply(fits[2:1], function(x) as.numeric(logLik(x)), numeric(1))))
  expect_equal(tab$AIC, unname(vapply(fits[2:1], AIC, numeric(1))))
  expect_equal(tab$BIC, unname(vapply(fits[2:1], BIC, numeric(1))))
  expect_false("L.Ratio" %in% names(tab))
  expect_null(attr(tab, "integration"))
  invisible(anova(fits$reduced, fits$full, nMC = 4096))
  expect_identical(serialize(fits, NULL), original)
})

test_that("multiple comparisons include survival covariates and pfs FALSE fits", {
  fits <- comparison_fixtures
  expect_null(fits$extended$Hessian)
  expect_true(is.list(fits$extended$dmats))
  tab <- anova(fits$reduced, fits$full, extended = fits$extended, nMC = 4096)
  expect_equal(nrow(tab), 3)
  expect_equal(tab$Df, c(NA, 1, 2))
  expect_equal(rownames(tab)[3], "extended")
  expect_true(all(is.finite(tab[["Pr(>Chisq)"]][-1])))
})

test_that("likelihood re-evaluation agrees with the existing joint likelihood", {
  dat <- small_joint_data()
  dat$T[dat$id == 3] <- 2.9 # censored before the first event
  fits <- c(comparison_fixtures,
    list(fit_small_joint(formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time),
                         formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id)),
         fit_small_joint(data = list(dat, subset(dat, !(id %% 3 == 0 & time == 2))),
                         formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time),
                         formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id))))
  for (fit in fits) {
    nodes <- matrix(randtoolbox::halton(fit$finalnMC, dim = sum(fit$dims$r),
                                      normal = TRUE), nrow = fit$finalnMC)
    ll <- joineRML:::mjoint_comparison_loglik(fit, nodes)
    expect_equal(ll[2], as.numeric(logLik(fit)), tolerance = 1e-10)
    expect_true(all(is.finite(ll)))
  }
})

test_that("likelihood re-evaluation agrees with independent scalar quadrature", {
  fit <- comparison_fixtures$full
  theta <- coef(fit)
  dat <- fit$data[[1]]
  surv <- fit$survData
  event.times <- sort(unique(surv$T[surv$status == 1]))
  expected <- sum(vapply(seq_len(nrow(surv)), function(i) {
    d <- dat[dat$id == surv$id[i], ]
    mu <- as.vector(model.matrix(fit$formLongFixed[[1]], d) %*% theta$beta)
    v <- theta$D[1, 1]
    s <- theta$sigma2[1]
    m <- nrow(d)
    fy <- mvtnorm::dmvnorm(d$y1, mu, diag(s, m) + matrix(v, m, m), log = TRUE)
    post.var <- 1 / (1 / v + m / s)
    post.mean <- post.var * sum(d$y1 - mu) / s
    H <- sum(theta$haz[event.times <= surv$T[i]])
    hazard <- if (surv$status[i] == 1) theta$haz[match(surv$T[i], event.times)] else 1
    integrand <- function(b) {
      eta <- theta$gamma[1] * b
      exp(surv$status[i] * (log(hazard) + eta) - H * exp(eta)) *
        dnorm(b, post.mean, sqrt(post.var))
    }
    fy + log(integrate(integrand, -Inf, Inf, rel.tol = 1e-9)$value)
  }, numeric(1)))
  nodes <- matrix(randtoolbox::halton(32768, dim = 1, normal = TRUE), ncol = 1)
  evaluated <- joineRML:::mjoint_comparison_loglik(fit, nodes)[2]
  expect_lt(abs(evaluated - expected), 0.01)
})

test_that("nesting is checked from column spaces rather than parameter counts", {
  fits <- comparison_fixtures
  expect_error(anova(fits$full, fits$reduced), "nested")
  expect_error(anova(fits$full, fits$full), "add at least one")
  other <- fit_small_joint(formLongFixed = list(y1 = y1 ~ I(time^2)))
  expect_error(anova(fits$full, other), "nested")
  expect_s3_class(anova(fits$full, other, test = FALSE), "anova")
  # Different contrasts span the same factor effects.
  d <- data.frame(g = factor(rep(letters[1:3], each = 4)), x = 1:12)
  x0 <- model.matrix(~ g, d, contrasts.arg = list(g = "contr.sum"))
  x1 <- model.matrix(~ g + x, d, contrasts.arg = list(g = "contr.treatment"))
  expect_true(joineRML:::mjoint_nested_matrix(x0, x1))
  expect_false(joineRML:::mjoint_nested_matrix(x1, x0))
  expect_false(joineRML:::mjoint_nested_matrix(x0, cbind(x1, x1[, 2])))
})

test_that("different observations are rejected even when sample sizes agree", {
  fits <- comparison_fixtures
  changes <- list(
    function(f) { f$data[[1]]$id <- as.character(f$data[[1]]$id); f$data[[1]]$id[1] <- "other"; f },
    function(f) { f$data[[1]]$time[1] <- 0.5; f },
    function(f) { f$dmats$l$yi[[1]][1] <- f$dmats$l$yi[[1]][1] + 1; f },
    function(f) { f$dmats$t$survdat2$T[1] <- f$dmats$t$survdat2$T[1] + 1; f },
    function(f) { f$dmats$t$survdat2$delta[1] <- 1 - f$dmats$t$survdat2$delta[1]; f },
    function(f) { f$formLongFixed[[1]][[2]] <- quote(log(y1)); f }
  )
  for (change in changes) {
    changed <- change(fits$full)
    for (test in c(FALSE, TRUE)) {
      expect_error(anova(fits$reduced, changed, test = test), "same subjects, responses")
    }
  }
  extra <- fits$full
  extra$data[[1]]$unused <- 1
  expect_s3_class(anova(fits$reduced, extra, test = FALSE), "anova")
})

test_that("non-nested survival covariates are rejected", {
  fit <- comparison_fixtures$extended
  other <- suppressMessages(mjoint(fit$formLongFixed, fit$formLongRandom,
    Surv(T, status) ~ group, data = fit$data, timeVar = "time",
    inits = list(gamma = c(0.1, 0.1)), control = fit$control))
  expect_error(anova(fit, other), "nested")
  expect_s3_class(anova(fit, other, test = FALSE), "anova")
})

test_that("random-effects changes do not receive an ordinary chi-squared test", {
  dat <- small_joint_data()
  dat$y1 <- dat$y1 + 0.8 * sin(dat$id) * dat$time
  intercept <- fit_small_joint(data = dat)
  slope <- fit_small_joint(data = dat, formLongRandom = list(y1 = ~ time | id))
  expect_error(anova(intercept, slope), "same random-effects and association")
  expect_s3_class(anova(intercept, slope, test = FALSE), "anova")
  nodes <- matrix(randtoolbox::halton(128, dim = 2, normal = TRUE), ncol = 2)
  expect_equal(joineRML:::mjoint_comparison_loglik(slope, nodes)[2],
               as.numeric(logLik(slope)), tolerance = 1e-10)
})

test_that("negative and unstable likelihood-ratio statistics withhold p-values", {
  fits <- comparison_fixtures
  poor <- fits$full
  poor$coefficients$beta[1] <- poor$coefficients$beta[1] + 10
  # A poorly maximised full model must not yield an absolute-value LRT.
  expect_warning(tab <- anova(fits$reduced, poor, nMC = 4096, tol = 1e6), "negative")
  expect_lt(tab$L.Ratio[2], 0)
  expect_true(is.na(tab[["Pr(>Chisq)"]][2]))
  expect_warning(tab <- anova(fits$reduced, fits$full, nMC = 128, tol = 1e-12), "unstable")
  expect_gt(tab$LR.change[2], 1e-12)
  expect_true(is.na(tab[["Pr(>Chisq)"]][2]))
})

test_that("invalid or unsuitable fitted objects give informative errors", {
  fits <- comparison_fixtures
  expect_error(anova(fits$full), "at least two")
  expect_error(anova(fits$full, 1), "All objects")
  for (value in list(NA, 1, "TRUE", c(TRUE, FALSE))) {
    expect_error(anova(fits$reduced, fits$full, test = value), "'test'")
  }
  for (value in list(NA, Inf, 1, 2.5, "128", c(2, 4))) {
    expect_error(anova(fits$reduced, fits$full, nMC = value), "'nMC'")
  }
  for (value in list(NA, Inf, 0, -1, "0.1", c(1, 2))) {
    expect_error(anova(fits$reduced, fits$full, tol = value), "'tol'")
  }
  bad <- fits$full
  bad$conv <- FALSE
  expect_error(anova(fits$reduced, bad), "must have converged")
  bad <- fits$full
  bad$dmats <- NULL
  expect_error(anova(fits$reduced, bad), "refit it")
  bad <- fits$full
  bad$log.lik <- NA_real_
  expect_error(anova(fits$reduced, bad), "finite log-likelihoods")
})

test_that("deterministic integration preserves R's random-number state", {
  fits <- comparison_fixtures
  set.seed(328)
  previous <- .Random.seed
  a <- anova(fits$reduced, fits$full, nMC = 4096)
  expect_identical(.Random.seed, previous)
  runif(3)
  b <- anova(fits$reduced, fits$full, nMC = 4096)
  expect_identical(a, b)
  previous <- .Random.seed
  on.exit(assign(".Random.seed", previous, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(anova(fits$reduced, fits$full, nMC = 4096))
  exists.after <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expect_false(exists.after)
})
