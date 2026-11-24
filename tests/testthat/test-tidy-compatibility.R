test_that("tidying identifies outcomes and handles transformations explicitly", {
  fit <- fit_small_joint(
    formLongFixed = list(first = y1 ~ time, second = y2 ~ time),
    formLongRandom = list(first = ~ 1 | id, second = ~ 1 | id))
  td <- tidy(fit, component = "all", conf.int = TRUE)
  expect_equal(td$component, c(rep("longitudinal", 4), rep("survival", 2)))
  expect_equal(td$outcome, c("first", "first", "second", "second", "first", "second"))
  expect_equal(td$effect, c(rep("fixed", 4), rep("association", 2)))
  expect_equal(td$estimate, unname(c(fit$coefficients$beta, fit$coefficients$gamma)))
  raw <- tidy(fit, conf.int = TRUE)
  hr <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  for (column in c("estimate", "conf.low", "conf.high")) {
    expect_equal(hr[[column]], exp(raw[[column]]))
  }
  expect_equal(hr$std.error, raw$std.error)
  expect_equal(hr$p.value, raw$p.value)
  expect_error(tidy(fit, component = "longitudinal", exponentiate = TRUE), "only supported")
  expect_error(tidy(fit, component = "all", exponentiate = TRUE), "only supported")
  expect_error(tidy(fit, ci = TRUE), "Unused argument.*ci")
  expect_error(tidy(fit, conf.level = 1), "strictly between")
  expect_error(tidy(fit, conf.level = NA), "strictly between")
  expect_error(tidy(fit, conf.int = NA), "TRUE or FALSE")
  expect_error(augment(fit, newdata = fit$data), "Unused argument.*newdata")
  expect_error(glance(fit, extra = TRUE), "Unused argument")
  bs <- structure(list(beta.se = rep(0.2, 4), gamma.se = rep(0.3, 2)), class = "bootSE")
  expect_equal(tidy(fit, bootSE = bs), tidy(fit, boot_se = bs))
  expect_equal(tidy(fit, boot_se = bs)$std.error, rep(0.3, 2))
  expect_error(tidy(fit, bootSE = bs, boot_se = bs), "only one")
  expect_error(tidy(fit, boot_se = list()), "not of class")
  unnamed <- fit
  names(unnamed$formLongFixed) <- NULL
  expect_equal(tidy(unnamed)$outcome, c("y1", "y2"))
  expect_true(all(c(".fitted_y1_0", ".fitted_y2_1") %in% names(augment(unnamed))))
})

test_that("augment matches supplied records without confusing fitted observations", {
  dat <- small_joint_data()
  fit <- fit_small_joint(dat)
  augmented <- augment(fit, data = dat)
  fitted.order <- match(paste(dat$id, dat$time),
                        paste(fit$data[[1]]$id, fit$data[[1]]$time))
  expect_equal(augmented$.fitted_y1_0, fitted(fit)[[1]][fitted.order], ignore_attr = TRUE)
  expect_equal(augmented$y1 - augmented$.fitted_y1_0, augmented$.resid_y1_0,
               ignore_attr = TRUE)
  shuffled <- dat[rev(seq_len(nrow(dat))), ]
  expect_equal(augment(fit, data = list(shuffled))$.fitted_y1_0,
               rev(augmented$.fitted_y1_0))
  shuffled$y1[1] <- 0
  expect_error(augment(fit, data = shuffled), "unchanged")
})

test_that("broom loading order preserves joineRML methods and arguments", {
  skip_if_not_installed("broom")
  fit <- fit_small_joint()
  fixture <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".log")
  on.exit(unlink(c(fixture, script, output)), add = TRUE)
  saveRDS(fit, fixture)
  for (packages in list(c("joineRML", "broom"), c("broom", "joineRML"))) {
    writeLines(c(
      paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
      sprintf("library(%s)", packages),
      paste0("fit <- readRDS(", deparse(fixture), ")"),
      "for (g in c('tidy', 'glance', 'augment')) {",
      "  stopifnot(identical(environment(getS3method(g, 'mjoint')), asNamespace('joineRML')))",
      "}",
      "bs <- structure(list(beta.se = rep(0.2, 2), gamma.se = 0.3), class = 'bootSE')",
      "stopifnot(identical(broom::tidy(fit, boot_se = bs), joineRML::tidy(fit, bootSE = bs)))",
      "stopifnot(identical(broom::glance(fit)$nobs, nobs(fit)))",
      "stopifnot(nrow(broom::augment(fit)) == 90L)",
      "stopifnot(nrow(broom::tidy(fit, component = 'all')) == 3L)"
    ), script)
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                      stdout = output, stderr = output)
    expect_equal(status, 0, info = paste(readLines(output), collapse = "\n"))
  }
})

test_that("optional table packages retain joint estimates and sample sizes", {
  skip_if_not_installed("modelsummary")
  skip_if_not_installed("gtsummary")
  config <- tempfile()
  old.config <- Sys.getenv("R_USER_CONFIG_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old.config)) Sys.unsetenv("R_USER_CONFIG_DIR") else
      Sys.setenv(R_USER_CONFIG_DIR = old.config)
    unlink(config, recursive = TRUE)
  }, add = TRUE)
  Sys.setenv(R_USER_CONFIG_DIR = config)
  fit <- fit_small_joint()
  estimates <- tidy(fit, conf.int = TRUE)
  table.model <- structure(list(tidy = estimates, glance = glance(fit)),
                           class = "modelsummary_list")
  expect_equal(modelsummary::get_estimates(table.model)$estimate, estimates$estimate)
  expect_equal(modelsummary::get_gof(table.model)$nobs, nobs(fit))
  table <- modelsummary::modelsummary(table.model, output = "data.frame",
    gof_map = data.frame(raw = c("nobs", "nevent", "AIC", "BIC"),
                        clean = c("Subjects", "Events", "AIC", "BIC"),
                        fmt = c(0, 0, 1, 1)))
  expect_s3_class(table, "data.frame")
  expect_equal(table[[4]][match(c("Subjects", "Events"), table$term)], c("30", "20"))
  combined <- tidy(fit, component = "all", conf.int = TRUE)
  formatted <- gtsummary::as_gtsummary(combined)
  expect_equal(formatted$table_body$estimate, combined$estimate)
  expect_equal(formatted$table_body$outcome, combined$outcome)
})
