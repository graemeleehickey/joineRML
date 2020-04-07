# test tidy, augment, glance methods for mjoint object (joineRML package)
context("Issue #59")
library(joineRML)

test_that("Example that was broken in #59 is not broken anymore", {
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  hvd <- hvd[, c(2:ncol(hvd), 1)]
  fit1 <- mjoint(
    formLongFixed = log.lvmi ~ time + age,
    formLongRandom = ~ time | num,
    formSurv = Surv(fuyrs, status) ~ age,
    data = hvd,
    timeVar = "time",
    control = list(nMCscale = 2, burnin = 5)
  )
  expect_s3_class(object = fit1, class = "mjoint")
})
