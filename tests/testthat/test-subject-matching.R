test_that("subject matching is independent of row order and factor levels", {
  dat <- small_joint_data()
  dat$x <- as.numeric(dat$id) / 30
  surv <- dat[!duplicated(dat$id), c("id", "T", "status", "x")]
  fit <- function(datasets, survival, inits = NULL) {
    suppressMessages(mjoint(
      formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time),
      formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id),
      formSurv = Surv(T, status) ~ x,
      data = datasets, survData = survival, timeVar = "time", inits = inits,
      control = list(type = "halton", nMC = 128, nMCmax = 128,
                     burnin = 1, mcmaxIter = 4, convCrit = "abs", tol0 = Inf)
    ))
  }
  reference <- fit(list(dat, dat), surv)
  set.seed(120)
  d1 <- dat[sample(nrow(dat)), ]
  d2 <- dat[sample(nrow(dat)), ]
  sd <- surv[sample(nrow(surv)), ]
  d1$id <- factor(d1$id, levels = rev(unique(dat$id)))
  d2$id <- factor(d2$id, levels = sample(unique(dat$id)))
  sd$id <- as.character(sd$id)
  reordered <- fit(list(d1, d2), sd)

  expect_equal(coef(reordered), coef(reference), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(reordered)), as.numeric(logLik(reference)),
               tolerance = 1e-8)
  expect_equal(reordered$Eb, reference$Eb, tolerance = 1e-8)
  expect_equal(reordered$Vb, reference$Vb, tolerance = 1e-8)
  expect_equal(vcov(reordered), vcov(reference), tolerance = 1e-8)

  # Repeat with unequal visit schedules and independently ordered ID factors.
  d1 <- subset(dat, !(id %% 3 == 0 & time == 2))
  d2 <- subset(dat, !(id %% 2 == 0 & time == 1))
  reference <- fit(list(d1, d2), surv, inits = list(gamma = c(0, 0.1, 0.1)))
  d1$id <- as.character(d1$id)
  d2$id <- factor(d2$id, levels = rev(unique(dat$id)))
  sd$id <- factor(sd$id, levels = sample(unique(dat$id)))
  reordered <- fit(list(d1[sample(nrow(d1)), ], d2[sample(nrow(d2)), ]), sd,
                   inits = list(gamma = c(0, 0.1, 0.1)))
  expect_equal(coef(reordered), coef(reference), tolerance = 1e-8)
  expect_equal(reordered$log.lik, reference$log.lik, tolerance = 1e-8)
  expect_equal(reordered$Eb, reference$Eb, tolerance = 1e-8)

  ids <- names(reordered$dmats$l$yi)
  expect_identical(as.character(reordered$survData$id), ids)
  for (x in c(reordered$dmats$l[c("Xi", "Zi", "Zit", "nik", "Xtyi", "XtZi")],
              reordered$dmats$t[c("V", "survdat2.list")], reordered$dmats$z)) {
    expect_identical(names(x), ids)
  }
})

test_that("invalid subject sets and missing survival values are rejected", {
  dat <- small_joint_data()
  surv <- dat[!duplicated(dat$id), c("id", "T", "status")]
  fit <- function(d = dat, s = surv, formula = Surv(T, status) ~ 1) {
    mjoint(y1 ~ time, ~ 1 | id, formula, d, survData = s, timeVar = "time")
  }
  expect_error(fit(s = rbind(surv, surv[1, ])), "exactly one record per subject")
  expect_error(fit(s = surv[-1, ]), "same subjects")
  extra <- surv[1, ]; extra$id <- 99
  expect_error(fit(s = rbind(surv, extra)), "same subjects")
  wrong <- surv; wrong$id[1] <- 99
  expect_error(fit(s = wrong), "same subjects")
  missing <- dat; missing$id[1] <- NA
  expect_error(fit(d = missing), "non-missing.*longitudinal")
  missing <- surv; missing$id[1] <- NA
  expect_error(fit(s = missing), "non-missing.*survData")
  expect_error(fit(s = surv[, c("T", "status")]), "IDs must be present")
  for (variable in c("T", "status")) {
    missing <- surv; missing[1, variable] <- NA
    expect_error(fit(s = missing), "Missing values in the survival model")
  }
  missing <- surv; missing$x <- seq_len(nrow(surv)); missing$x[1] <- NA
  expect_error(fit(s = missing, formula = Surv(T, status) ~ x),
               "Missing values in the survival model")
})

test_that("censoring the first subject before the first event is supported", {
  dat <- small_joint_data()
  dat$T[dat$id == 1] <- 2.5
  dat$status[dat$id == 1] <- 0
  fit <- fit_small_joint(data = dat)
  expect_true(fit$conv)
  expect_equal(fit$survData$status[fit$survData$id == "1"], 0)
  expect_equal(fit$dmats$t$survdat2.list[["1"]]$tj.ind, 0)
  expect_true(is.finite(fit$log.lik))
  theta <- coef(fit)
  l <- fit$dmats$l
  m <- nrow(l$Xi[["1"]])
  # With no baseline hazard before censoring, this subject's posterior is the
  # Gaussian posterior conditional on the longitudinal observations alone.
  post.var <- 1 / (1 / theta$D[1, 1] + m / theta$sigma2[1])
  post.mean <- post.var * sum(l$yi[["1"]] - l$Xi[["1"]] %*% theta$beta) / theta$sigma2[1]
  nodes <- as.vector(randtoolbox::halton(fit$finalnMC, dim = 1, normal = TRUE))
  row <- match("1", names(l$yi))
  expect_equal(unname(fit$Eb[row, 1]),
               unname(post.mean + sqrt(post.var) * mean(nodes)), tolerance = 1e-10)
})
