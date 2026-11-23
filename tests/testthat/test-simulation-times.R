test_that("Gompertz times include the exponential limit", {
  simulate <- function(slope) {
    set.seed(45)
    suppressMessages(simData(n = 100, model = "intslope", theta0 = -1,
      theta1 = slope, gamma.x = c(0, 0), gamma.y = c(0, 0),
      censoring = FALSE, truncation = FALSE))$survdat$survtime
  }
  t0 <- simulate(0)
  expect_true(all(is.finite(t0) & t0 > 0))
  rate <- exp(-1)
  for (slope in c(-0.2, -1e-10, 1e-10, 0.2)) {
    times <- simulate(slope)
    finite <- is.finite(times)
    expect_false(anyNA(times))
    expect_true(all(times > 0))
    # H_slope(T_slope) must equal H_0(T_0) for the same uniform draws.
    expect_equal(rate * expm1(slope * times[finite]) / slope,
                 rate * t0[finite], tolerance = 1e-10)
    if (slope < 0) {
      expect_equal(is.infinite(times), t0 >= -1 / slope)
    } else {
      expect_true(all(finite))
    }
    if (abs(slope) < 1e-9) {
      expect_equal(times, t0, tolerance = 1e-8)
    }
  }
})

test_that("simulated censoring status reflects the generated event time", {
  simulate <- function(censoring = FALSE, truncation = FALSE, trunctime = 4.1) {
    set.seed(45)
    suppressMessages(simData(n = 100, model = "int", theta0 = -1,
      gamma.x = c(0, 0), gamma.y = c(0, 0), censlam = 0.4,
      censoring = censoring, truncation = truncation, trunctime = trunctime))
  }
  events <- simulate()$survdat$survtime
  for (random.censoring in c(FALSE, TRUE)) {
    uncapped <- simulate(censoring = random.censoring)$survdat
    expect_equal(uncapped$cens, as.integer(uncapped$survtime == events))
    cap <- min(events) / 2
    capped <- simulate(censoring = random.censoring, truncation = TRUE,
                       trunctime = cap)$survdat
    expect_equal(capped$survtime, pmin(uncapped$survtime, cap))
    expect_equal(capped$cens, rep(0L, length(events)))
  }
  # Ties at a finite administrative censoring time are observed events.
  tied <- simulate(truncation = TRUE, trunctime = events[1])$survdat
  expect_equal(tied$cens, as.integer(events <= events[1]))
})

test_that("infinite Gompertz event times are not labelled as observed events", {
  set.seed(45)
  simulated <- suppressMessages(simData(n = 100, model = "intslope",
    theta0 = -1, theta1 = -0.2, gamma.x = c(0, 0), gamma.y = c(0, 0),
    censoring = FALSE, truncation = FALSE))
  expect_true(any(is.infinite(simulated$survdat$survtime)))
  expect_equal(simulated$survdat$cens,
               as.integer(is.finite(simulated$survdat$survtime)))
})

test_that("no-event data are returned faithfully and cannot be fitted", {
  set.seed(11)
  simulated <- suppressMessages(simData(n = 30, model = "int", theta0 = -100,
    censoring = FALSE, truncation = TRUE, trunctime = 1))
  expect_equal(simulated$survdat$cens, rep(0L, 30))
  expect_equal(simulated$survdat$survtime, rep(1, 30))
  expect_error(mjoint(Y.1 ~ time, ~ 1 | id, Surv(survtime, cens) ~ 1,
    data = simulated$longdat, survData = simulated$survdat, timeVar = "time"),
    "at least one observed event")
})
