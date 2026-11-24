test_that("residual summaries agree in scale and identify outcomes", {
  for (K in 1:2) {
    fit <- fit_small_joint(
      formLongFixed = list(first = y1 ~ time, second = y2 ~ time)[seq_len(K)],
      formLongRandom = list(first = ~ 1 | id, second = ~ 1 | id)[seq_len(K)],
      pfs = FALSE)
    for (named in c(TRUE, FALSE)) {
      if (!named) names(fit$formLongFixed) <- NULL
      variances <- coef(fit)$sigma2
      expect_equal(unname(sigma(fit)^2), unname(variances))
      expect_identical(summary(fit)$sigma, sigma(fit))
      columns <- paste0("sigma_", seq_len(K))
      expect_equal(as.numeric(glance(fit)[1, columns]), unname(sigma(fit)))
      expect_false(any(grepl("^sigma2_", names(glance(fit)))))

      # Known non-unit variances distinguish SDs from variances in the output.
      fit$coefficients$sigma2 <- setNames(c(0.25, 4)[seq_len(K)],
                                         paste0("sigma2_", seq_len(K)))
      expected <- setNames(c(0.5, 2)[seq_len(K)],
        if (named) c("first", "second")[seq_len(K)] else columns)
      expect_identical(sigma(fit), expected)
      expect_identical(summary(fit)$sigma, expected)
      expect_equal(as.numeric(glance(fit)[1, columns]), unname(expected))
      printed <- capture.output(print(summary(fit)))
      expect_true(any(grepl("Residual standard deviations:", printed, fixed = TRUE)))
      expect_false(any(grepl("sigma2_|Residual standard errors", printed)))
      expect_identical(coef(fit)$sigma2,
                        setNames(c(0.25, 4)[seq_len(K)], paste0("sigma2_", seq_len(K))))
      fit$coefficients$sigma2 <- variances
    }
  }
})

test_that("bootstrap printing preserves residual variance estimates and uncertainty", {
  # Controlled bootstrap results make a change of scale detectable in each column.
  bs <- structure(list(
    coefficients = list(beta = c("(Intercept)_1" = 1), gamma = c(gamma_1 = 0.1),
                        sigma2 = c(sigma2_1 = 0.25, sigma2_2 = 4)),
    beta.se = 0.1, gamma.se = 0.2,
    beta.ci = matrix(c(0.8, 1.2), ncol = 1),
    gamma.ci = matrix(c(-0.3, 0.5), ncol = 1),
    sigma2.se = c(0.03, 0.4),
    sigma2.ci = matrix(c(0.16, 0.36, 3, 5), nrow = 2),
    ci = 0.95, boot.time = as.difftime(1, units = "secs"), conv = 3, nboot = 3
  ), class = "bootSE")
  printed <- capture.output(returned <- print(bs))
  expect_identical(returned, bs)
  expect_true(any(grepl("Residual variances", printed, fixed = TRUE)))
  expect_false(any(grepl("Residual SE", printed, fixed = TRUE)))
  for (k in 1:2) {
    row <- printed[grepl(paste0("sigma2_", k), printed, fixed = TRUE)]
    values <- as.numeric(utils::tail(strsplit(trimws(row), "\\s+")[[1]], 4))
    expect_equal(values, c(bs$coefficients$sigma2[k], bs$sigma2.se[k],
                           bs$sigma2.ci[, k]), ignore_attr = TRUE)
  }
})
