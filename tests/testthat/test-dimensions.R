test_that("fitted values and residuals retain observations and their order", {
  dat <- small_joint_data()
  for (K in 1:2) {
    for (equal.visits in c(TRUE, FALSE)) {
      datasets <- rep(list(dat), K)
      if (!equal.visits) {
        datasets[[1]] <- subset(dat, !(id %% 3 == 0 & time == 2))
        if (K == 2) {
          datasets[[2]] <- subset(dat, !(id %% 2 == 0 & time == 1))
        }
      }
      fit <- fit_small_joint(
        data = datasets,
        formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time)[seq_len(K)],
        formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id)[seq_len(K)]
      )
      for (level in 0:1) {
        fv <- fitted(fit, level = level)
        rv <- residuals(fit, level = level)
        expect_named(fv, names(fit$formLongFixed))
        expect_named(rv, names(fit$formLongFixed))
        for (k in seq_len(K)) {
          d <- fit$data[[k]]
          beta.k <- coef(fit)$beta[(2 * k - 1):(2 * k)]
          expected <- beta.k[1] + beta.k[2] * d$time
          if (level == 1) {
            subjects <- unique(fit$data[[1]]$id)
            expected <- expected + fit$Eb[match(d$id, subjects), k]
          }
          observed <- d[[paste0("y", k)]]
          expect_null(dim(fv[[k]]))
          expect_null(dim(rv[[k]]))
          expect_length(fv[[k]], nrow(d))
          expect_length(rv[[k]], nrow(d))
          expect_equal(as.numeric(fv[[k]]), as.numeric(expected))
          expect_equal(as.numeric(rv[[k]]), as.numeric(observed - expected))
          expect_equal(as.numeric(fv[[k]] + rv[[k]]), observed)
        }
      }
    }
  }
})

test_that("intercept-only fixed effects work through post-fit calculations", {
  for (pfs in c(FALSE, TRUE)) {
    fit <- fit_small_joint(formLongFixed = list(y1 = y1 ~ 1), pfs = pfs)
    expect_true(fit$conv)
    expect_length(coef(fit)$beta, 1)
    expect_equal(ncol(fit$history), 4)
    expect_true("long_(Intercept)_1" %in% rownames(fit$history))
    expect_true(is.finite(as.numeric(logLik(fit))))
    if (pfs) {
      expect_equal(dim(vcov(fit)), c(4L, 4L))
      expect_true(all(is.finite(vcov(fit))))
      expect_length(fitted(fit)[[1]], nrow(fit$data[[1]]))
    } else {
      expect_null(fit$Hessian)
    }
  }

  fit <- fit_small_joint(
    formLongFixed = list(y1 = y1 ~ 1, y2 = y2 ~ time),
    formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id)
  )
  expect_true(fit$conv)
  expect_equal(fit$dims$p, c(1L, 2L))
  expect_length(coef(fit)$beta, 3)
  expect_true(all(is.finite(vcov(fit))))
})

test_that("bootstrap summaries retain single-parameter dimensions", {
  for (intercept.only in c(FALSE, TRUE)) {
    fixed <- if (intercept.only) list(y1 = y1 ~ 1) else list(y1 = y1 ~ time)
    fit <- fit_small_joint(formLongFixed = fixed)
    for (nboot in c(1, 3)) {
      set.seed(129)
      bs <- bootSE(fit, nboot = nboot, progress = FALSE)
      expect_equal(bs$conv, nboot)
      expect_equal(dim(bs$beta.ci), c(2L, length(coef(fit)$beta)))
      expect_equal(dim(bs$gamma.ci), c(2L, 1L))
      expect_equal(dim(bs$sigma2.ci), c(2L, 1L))
      expect_named(bs$gamma.se, "gamma_1")
      expect_equal(colnames(bs$beta.ci), names(coef(fit)$beta))
      if (nboot == 1) {
        expect_true(all(is.na(c(bs$beta.se, bs$gamma.se, bs$sigma2.se))))
      } else {
        expect_true(all(is.finite(c(bs$beta.se, bs$gamma.se, bs$sigma2.se))))
      }
    }
  }
})

test_that("one distinct event time supports fitting and hazard histories", {
  dat <- small_joint_data()
  dat$T <- 4
  for (K in 1:2) {
    for (pfs in c(FALSE, TRUE)) {
      for (early.censoring in c(FALSE, TRUE)) {
        dat$T[dat$status == 0] <- if (early.censoring) 3 else 4
        fit <- fit_small_joint(dat, pfs = pfs,
          formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time)[seq_len(K)],
          formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id)[seq_len(K)])
        expect_true(fit$conv)
        expect_true(all(is.finite(unlist(coef(fit)))))
        expect_true(is.finite(as.numeric(logLik(fit))))
        expect_equal(baseHaz(fit, centered = FALSE)$time, 4)
        expect_length(coef(fit)$haz, 1)
        expect_gt(coef(fit)$haz, 0)
        expect_equal(as.numeric(tail(fit$history["haz_1", ], 1)), coef(fit)$haz)
        if (pfs) {
          expect_equal(dim(fit$Eb), c(30L, K))
          expect_true(all(is.finite(fit$Hessian)))
        }
      }
    }
  }
})

test_that("extraction explains unavailable posterior random effects", {
  for (pfs in c(FALSE, TRUE)) {
    for (iterations in c(2, 4)) {
      fit <- fit_small_joint(pfs = pfs, mcmaxIter = iterations)
      dat <- fit$data[[1]]
      expected <- as.vector(model.matrix(fit) %*%
                              coef(fit, process = "Longitudinal", k = 1))
      expect_equal(as.numeric(fitted(fit, level = 0)[[1]]), expected)
      expect_equal(as.numeric(residuals(fit, level = 0)[[1]]), dat$y1 - expected)
      if (!pfs || !fit$conv) {
        message <- "Posterior random effects.*pfs = TRUE.*convergence"
        expect_error(fitted(fit, level = 1), message)
        expect_error(residuals(fit, level = 1), message)
        expect_error(augment(fit), message)
        expect_error(ranef(fit), message)
        expect_error(ranef(fit, postVar = TRUE), message)
      } else {
        expect_length(fitted(fit, level = 1)[[1]], nrow(dat))
        expect_length(residuals(fit, level = 1)[[1]], nrow(dat))
        expect_equal(nrow(augment(fit)), nrow(dat))
        expect_equal(nrow(ranef(fit)), fit$dims$n)
        expect_equal(attr(ranef(fit, postVar = TRUE), "postVar"), fit$Vb)
      }
    }
  }
})
