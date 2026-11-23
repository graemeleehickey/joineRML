# Short fits exercise the full sampling and post-fit calculations; the loose
# convergence settings here are only for testing reproducibility.
seed_fit_args <- function(type = "antithetic", pfs = TRUE, K = 1) {
  list(
    formLongFixed = list(y1 = y1 ~ time, y2 = y2 ~ time)[seq_len(K)],
    formLongRandom = list(y1 = ~ 1 | id, y2 = ~ 1 | id)[seq_len(K)],
    formSurv = Surv(T, status) ~ 1,
    data = small_joint_data(),
    timeVar = "time",
    inits = list(gamma = rep(0.1, K)),
    pfs = pfs,
    control = list(type = type, nMC = 128, nMCmax = 128,
                   burnin = 1, mcmaxIter = 4, convCrit = "abs",
                   tol0 = Inf, tol.em = 1e-3)
  )
}

seed_fit_results <- function(fit) {
  fit[c("coefficients", "history", "log.lik", "ll.hx", "Hessian", "Eb", "Vb")]
}

test_that("seed controls both R and C++ sampling throughout the fit", {
  for (type in c("montecarlo", "antithetic")) {
    for (K in 1:2) {
      for (pfs in c(FALSE, TRUE)) {
        args <- seed_fit_args(type, pfs, K)
        set.seed(729)
        previous <- .Random.seed
        fit1 <- suppressMessages(do.call(mjoint, c(args, list(seed = 153))))
        expect_identical(.Random.seed, previous)
        expect_true(fit1$conv)
        expect_identical(fit1$seed, 153L)
        expect_identical(fit1$call$seed, 153L)
        if (pfs) {
          expect_true(all(is.finite(fit1$Hessian)))
        }

        runif(5)
        previous <- .Random.seed
        fit2 <- suppressMessages(do.call(mjoint, c(args, list(seed = 153L))))
        expect_identical(.Random.seed, previous)
        expect_identical(seed_fit_results(fit1), seed_fit_results(fit2))

        fit3 <- suppressMessages(do.call(mjoint, c(args, list(seed = 154))))
        expect_false(identical(fit1$history, fit3$history))

        set.seed(153)
        external <- suppressMessages(do.call(mjoint, args))
        expect_identical(seed_fit_results(fit1), seed_fit_results(external))
      }
    }
  }
})

test_that("omitted and NULL seeds retain the existing random-number sequence", {
  args <- seed_fit_args()
  set.seed(17)
  previous <- .Random.seed
  omitted <- suppressMessages(do.call(mjoint, args))
  following <- .Random.seed
  expect_false(identical(previous, following))
  expect_null(omitted$seed)
  set.seed(17)
  explicit <- suppressMessages(do.call(mjoint, c(args, list(seed = NULL))))
  expect_identical(seed_fit_results(omitted), seed_fit_results(explicit))
  expect_identical(.Random.seed, following)
})

test_that("quasi-Monte Carlo fits are reproducible with a supplied seed", {
  for (type in c("sobol", "halton")) {
    args <- seed_fit_args(type)
    # randtoolbox versions with scrambling disabled warn on every Sobol call.
    fit1 <- suppressWarnings(suppressMessages(
      do.call(mjoint, c(args, list(seed = 21)))))
    runif(7)
    fit2 <- suppressWarnings(suppressMessages(
      do.call(mjoint, c(args, list(seed = 21)))))
    expect_true(fit1$conv)
    expect_identical(seed_fit_results(fit1), seed_fit_results(fit2))
  }
})

test_that("the seed is retained when updating a fitted model", {
  args <- seed_fit_args()
  chosen.seed <- 37
  fit <- suppressMessages(do.call(mjoint, c(args, list(seed = chosen.seed))))
  chosen.seed <- 38
  repeated <- suppressMessages(update(fit))
  expect_identical(seed_fit_results(fit), seed_fit_results(repeated))
  expect_identical(repeated$seed, 37L)
})

test_that("invalid seeds fail without changing the random-number state", {
  set.seed(52)
  previous <- .Random.seed
  for (seed in list(NA, NA_real_, Inf, NaN, -1, 1.5, "1", TRUE,
                    numeric(), c(1, 2), 1 + 1i, .Machine$integer.max + 1)) {
    expect_error(mjoint(seed = seed), "'seed' must be NULL or a single")
    expect_identical(.Random.seed, previous)
  }
})

test_that("the random-number state is restored after an error", {
  args <- seed_fit_args()
  args$data$id[1] <- NA
  set.seed(813)
  previous <- .Random.seed
  expect_error(do.call(mjoint, c(args, list(seed = 19))), "Subject IDs")
  expect_identical(.Random.seed, previous)
})

test_that("a supplied seed respects an alternative R random-number generator", {
  args <- seed_fit_args()
  previous.kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(previous.kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(381)
  previous <- .Random.seed
  fit <- suppressMessages(do.call(mjoint, c(args, list(seed = 123))))
  expect_identical(.Random.seed, previous)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  set.seed(123)
  external <- suppressMessages(do.call(mjoint, args))
  expect_identical(seed_fit_results(fit), seed_fit_results(external))
})

test_that("a supplied seed works when no previous random-number state exists", {
  args <- seed_fit_args()
  previous <- .Random.seed
  on.exit(assign(".Random.seed", previous, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  fit <- suppressMessages(do.call(mjoint, c(args, list(seed = 0))))
  exists.after.fit <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  args$data$id[1] <- NA
  error <- tryCatch(do.call(mjoint, c(args, list(seed = .Machine$integer.max))),
                     error = identity)
  exists.after.error <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  # Capture the state before testthat creates its own random identifiers.
  expect_true(fit$conv)
  expect_false(exists.after.fit)
  expect_match(conditionMessage(error), "Subject IDs")
  expect_false(exists.after.error)
})
