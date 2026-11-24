prediction_basis_data <- function() {
  set.seed(2314)
  n <- 30
  visits <- seq(0, 2, length.out = 6)
  id <- rep(seq_len(n), each = length(visits))
  time <- rep(visits, n)
  data.frame(id = id, time = time,
    y1 = 2 + 0.3 * time + rnorm(n)[id] + rnorm(n, sd = 0.4)[id] * time +
      rnorm(n, sd = 0.3)[id] * time^2 + rnorm(length(id), sd = 0.25),
    T = rep(3 + rexp(n), each = length(visits)),
    status = rep(rep(c(1, 1, 0), length.out = n), each = length(visits)))
}

test_that("fitting and prediction retain the same random-effects basis", {
  dat <- prediction_basis_data()
  cases <- list(
    list(formula = ~ poly(time, 2) | id, basis = poly(dat$time, 2)),
    list(formula = ~ stats::poly(time, 1) | id, basis = poly(dat$time, 1)),
    list(formula = ~ splines::ns(time, df = 2) | id,
         basis = splines::ns(dat$time, df = 2)),
    list(formula = ~ splines::bs(time, degree = 2, df = 2,
                                Boundary.knots = c(0, 10)) | id,
         basis = splines::bs(dat$time, degree = 2, df = 2,
                             Boundary.knots = c(0, 10))),
    list(formula = ~ scale(time) | id, basis = scale(dat$time)))
  for (case in cases) {
    fit <- fit_small_joint(dat, formLongRandom = list(y1 = case$formula), pfs = FALSE)
    expected <- function(time) {
      if (inherits(case$basis, c("poly", "ns", "bs"))) {
        cbind(1, predict(case$basis, time))
      } else {
        cbind(1, (time - attr(case$basis, "scaled:center")) /
                    attr(case$basis, "scaled:scale"))
      }
    }
    expect_equal(unname(do.call(rbind, fit$dmats$l$Zik.list[[1]])),
                  unname(expected(fit$data[[1]]$time)), ignore_attr = TRUE)
    for (i in c(1, 4, 13)) {
      times <- fit$dmats$t$tj[seq_len(max(fit$dmats$t$survdat2$tj.ind[i], 1))]
      expect_equal(unname(fit$dmats$z$Zi.fail[[i]]), unname(expected(times)))
    }
    nd <- subset(fit$data[[1]], id == 1)
    history <- joineRML:::process_newdata(fit, nd)
    expect_equal(unname(history$Z), unname(expected(nd$time)))
    future <- nd[1, , drop = FALSE]
    future$time <- 3
    one <- joineRML:::process_newdata(fit, future)
    many <- future[rep(1, 3), ]
    many$time <- c(2.5, 3, 3.5)
    batch <- joineRML:::process_newdata(fit, many)
    expect_equal(unname(one$Z), unname(expected(3)))
    expect_equal(unname(one$Z), unname(batch$Z[2, , drop = FALSE]))
    expect_equal(unname(one$Z.fail),
                  unname(expected(fit$dmats$t$tj[seq_len(max(one$tj.ind, 1))])))
    dl <- dynLong(fit, nd, u = 3, progress = FALSE)
    ds <- dynSurv(fit, nd, u = 3.5, progress = FALSE)
    expect_true(all(is.finite(as.matrix(dl$pred[[1]]))))
    expect_true(all(is.finite(as.matrix(ds$pred))))
  }
})

test_that("dynamic predictions use fitted factor contrasts in both sub-models", {
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old), add = TRUE)
  for (coding in c("contr.treatment", "contr.sum", "custom")) {
    dat <- small_joint_data()
    dat$group <- factor(letters[1 + ((dat$id - 1) %/% 2 %% 3)])
    dat$age <- dat$id / 30
    if (coding == "custom") {
      contrasts(dat$group) <- matrix(c(1, 0, -1, 0, 1, -1), 3,
                                     dimnames = list(letters[1:3], c("a", "b")))
    } else {
      options(contrasts = c(coding, "contr.poly"))
    }
    fit <- fit_small_joint(dat, formLongFixed = list(y1 = y1 ~ time * group),
                          formSurv = Surv(T, status) ~ group + age)
    nd <- droplevels(subset(fit$data[[1]], id == 1))
    # Prediction does not require an event response in the separate data.
    sd <- nd[1, c("group", "age"), drop = FALSE]
    history <- joineRML:::process_newdata(fit, nd, newSurvData = sd)
    prediction <- function(type) {
      set.seed(982)
      list(long = dynLong(fit, nd, newSurvData = sd, u = 3, type = type,
                          M = 20, progress = FALSE)$pred,
           surv = dynSurv(fit, nd, newSurvData = sd, u = 3.5, type = type,
                          M = 20, progress = FALSE)$pred)
    }
    first <- prediction("first-order")
    simulated <- prediction("simulated")
    options(contrasts = c("contr.helmert", "contr.SAS"))
    expect_identical(joineRML:::process_newdata(fit, nd, newSurvData = sd), history)
    expect_identical(prediction("first-order"), first)
    expect_identical(prediction("simulated"), simulated)
    older <- fit
    older$designs <- NULL
    older$model.frames <- NULL
    expect_equal(dynLong(older, nd, newSurvData = sd, u = 3)$pred, first$long)
    expect_equal(dynSurv(older, nd, newSurvData = sd, u = 3.5)$pred, first$surv)
    expected <- fit$sfit$x[1, ] - fit$dmats$t$xcenter
    expect_equal(unname(history$v), unname(expected))
    # Restoring levels must also handle a different ordering in the new data.
    reordered <- nd
    reordered$group <- factor(as.character(nd$group), levels = rev(levels(dat$group)))
    expect_identical(joineRML:::process_newdata(fit, reordered)$X, history$X)
    nd$group <- factor("new")
    expect_error(dynLong(fit, nd, u = 3), "new level")
    sd$group <- factor("new")
    expect_error(dynSurv(fit, subset(fit$data[[1]], id == 1),
                         newSurvData = sd, u = 3.5), "new level")
  }
})

test_that("older fits require agreement with their fitted random-effects matrices", {
  dat <- prediction_basis_data()
  fit <- fit_small_joint(dat, formLongRandom = list(y1 = ~ time | id), pfs = FALSE)
  nd <- subset(fit$data[[1]], id == 1)
  older <- fit
  older$designs <- NULL
  older$model.frames <- NULL
  expect_equal(dynLong(older, nd, u = 3)$pred, dynLong(fit, nd, u = 3)$pred)
  expect_equal(dynSurv(older, nd, u = 3.5)$pred, dynSurv(fit, nd, u = 3.5)$pred)
  older$dmats <- NULL
  expect_error(dynLong(older, nd, u = 3), "basis cannot be verified.*refit")

  fit <- fit_small_joint(dat, formLongRandom = list(y1 = ~ poly(time, 1) | id), pfs = FALSE)
  older <- fit
  older$designs <- NULL
  # Reproduce the old calculation: a new polynomial basis on the event times.
  event.times <- lapply(fit$dmats$t$survdat2$tj.ind, function(n) {
    fit$dmats$t$tj[seq_len(max(n, 1))]
  })
  old.basis <- model.matrix(~ poly(time, 1), data.frame(time = unlist(event.times)))
  ends <- cumsum(lengths(event.times))
  starts <- c(1L, head(ends, -1) + 1L)
  older$dmats$z$Zi.fail <- lapply(seq_along(ends), function(i) {
    old.basis[starts[i]:ends[i], , drop = FALSE]
  })
  expect_error(dynLong(older, nd, u = 3), "inconsistent random-effects bases.*refit")
  expect_error(dynSurv(older, nd, u = 3.5), "inconsistent random-effects bases.*refit")
})

test_that("prediction rejects changed design columns", {
  fit <- fit_small_joint()
  fit$designs$fixed[[1]]$columns <- rev(fit$designs$fixed[[1]]$columns)
  expect_error(dynLong(fit, subset(fit$data[[1]], id == 1), u = 3),
                "prediction design columns do not match")
})

test_that("outcomes with different schedules retain separate fitted bases", {
  dat <- prediction_basis_data()
  dat2 <- subset(dat, !(id %% 2 == 0 & time == 0.4))
  names(dat2)[names(dat2) == "time"] <- "visit"
  names(dat2)[names(dat2) == "y1"] <- "y2"
  dat2$y2 <- 0.5 * dat2$y2 + sin(dat2$id)
  fit <- suppressMessages(mjoint(
    list(y1 = y1 ~ time, y2 = y2 ~ visit),
    list(y1 = ~ poly(time, 1) | id, y2 = ~ scale(visit) | id),
    Surv(T, status) ~ 1, data = list(dat, dat2), timeVar = c("time", "visit"),
    inits = list(gamma = c(0.1, 0.1)), pfs = FALSE,
    control = list(type = "halton", nMC = 128, nMCmax = 128, burnin = 1,
                   mcmaxIter = 4, convCrit = "abs", tol0 = Inf, tol.em = 1e-3)))
  nd <- lapply(fit$data, function(d) subset(d, id == 2))
  basis1 <- poly(dat$time, 1)
  basis2 <- scale(dat2$visit)
  history <- joineRML:::process_newdata(fit, nd)
  expect_equal(unname(history$Zk[[1]]), unname(cbind(1, predict(basis1, nd[[1]]$time))),
                ignore_attr = TRUE)
  expect_equal(unname(history$Zk[[2]]),
                unname(cbind(1, (nd[[2]]$visit - attr(basis2, "scaled:center")) /
                              attr(basis2, "scaled:scale"))), ignore_attr = TRUE)
  i <- match("2", names(fit$dmats$z$Zi.fail))
  times <- fit$dmats$t$tj[seq_len(max(fit$dmats$t$survdat2$tj.ind[i], 1))]
  expected <- as.matrix(Matrix::bdiag(cbind(1, predict(basis1, times)),
    cbind(1, (times - attr(basis2, "scaled:center")) / attr(basis2, "scaled:scale"))))
  expect_equal(fit$dmats$z$Zi.fail[[i]], expected)
  long <- dynLong(fit, nd, u = 3)
  expect_length(long$pred, 2)
  expect_true(all(is.finite(unlist(long$pred))))
  expect_true(all(is.finite(as.matrix(dynSurv(fit, nd, u = 3.5)$pred))))
})
