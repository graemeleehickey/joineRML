reprex::reprex({
  library(joineRML)
  data(heart.valve)
  hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
  # Put ID last
  hvd <- hvd[, c(2:ncol(hvd), 1)]

  # Fit model
  set.seed(1)
  fit1 <- mjoint(
    formLongFixed = log.lvmi ~ time + age,
    formLongRandom = ~ time | num,
    formSurv = Surv(fuyrs, status) ~ age,
    data = hvd,
    timeVar = "time",
    control = list(nMCscale = 2, burnin = 5)
  )
})
