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
