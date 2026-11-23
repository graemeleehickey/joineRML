# A local copy lets tests specify a known target without changing package code.
metropolis_for_target <- function(log.target, uniform = NULL) {
  kernel <- joineRML:::b_metropolis
  env <- new.env(parent = environment(kernel))
  env$logpb <- function(b, theta, data) log.target(b)
  if (!is.null(uniform)) {
    env$runif <- function(n) uniform
  }
  environment(kernel) <- env
  kernel
}

test_that("identical target and proposal densities accept every proposal", {
  kernel <- metropolis_for_target(function(b) stats::dt(b, df = 4, log = TRUE))
  set.seed(184)
  accepted <- replicate(200, kernel(NULL, 0, matrix(1), 0, NULL)$accept)
  expect_true(all(accepted == 1))
})

test_that("Metropolis decisions agree with the analytical density ratio", {
  target <- function(b) stats::dnorm(b, mean = 0.3, sd = 0.7, log = TRUE)
  delta <- 1.1
  variance <- 0.7
  proposal <- function(b) stats::dt((b - delta) / sqrt(variance), df = 4,
                                    log = TRUE) - log(sqrt(variance))
  for (current in c(-2, -0.5, 1, 4)) {
    for (seed in 1:4) {
      set.seed(seed)
      candidate <- as.vector(mvtnorm::rmvt(1, delta = delta,
                                          sigma = matrix(variance), df = 4))
      alpha <- min(1, exp(target(candidate) - target(current) +
                          proposal(current) - proposal(candidate)))
      for (u in c(alpha / 2, (alpha + 1) / 2)) {
        kernel <- metropolis_for_target(target, uniform = u)
        set.seed(seed)
        result <- kernel(NULL, delta, matrix(variance), current, NULL)
        accept <- u <= alpha
        expect_equal(result$accept, as.numeric(accept))
        expect_equal(result$b.curr, if (accept) candidate else current)
      }
    }
  }
})

test_that("random-intercept posterior moments agree with numerical integration", {
  theta <- list(D = matrix(0.8), beta = c(0.2, 0.1), sigma2 = 0.4,
                gamma = 0.7, haz = c(0.05, 0.12))
  d <- list(r = 1, q = 0, K = 1, nk = 3, y = c(0.5, 1.3, 0.9),
            X = cbind(1, 0:2), Z = matrix(1, 3, 1), v = NULL,
            Z.fail = matrix(1, 2, 1), IW.fail = diag(2), tj.ind = 2)
  log.target <- function(b) {
    dnorm(b, 0, sqrt(0.8), log = TRUE) +
      sum(dnorm(d$y, 0.2 + 0.1 * (0:2) + b, sqrt(0.4), log = TRUE)) -
      0.17 * exp(0.7 * b)
  }
  mode <- optimize(log.target, c(-5, 5), maximum = TRUE)$maximum
  density <- function(b) vapply(b, function(x) exp(log.target(x)), numeric(1))
  mass <- integrate(density, -Inf, Inf, rel.tol = 1e-10)$value
  moments <- vapply(1:2, function(power) {
    integrate(function(b) b^power * density(b), -Inf, Inf, rel.tol = 1e-10)$value / mass
  }, numeric(1))
  sigma <- matrix(2 / (1 / 0.8 + 3 / 0.4 + 0.17 * 0.7^2 * exp(0.7 * mode)))
  set.seed(2026)
  current <- mode
  draws <- numeric(7000)
  for (i in seq_along(draws)) {
    current <- joineRML:::b_metropolis(theta, mode, sigma, current, d)$b.curr
    draws[i] <- current
  }
  draws <- draws[-seq_len(1000)]
  for (power in 1:2) {
    batches <- colMeans(matrix(draws^power, nrow = 100))
    mcse <- sd(batches) / sqrt(length(batches))
    expect_lt(abs(mean(draws^power) - moments[power]), 6 * mcse + 0.002)
  }
})
