# Dynamic predictions for the time-to-event data sub-model

Calculates the conditional time-to-event distribution for a *new*
subject from the last observation time given their longitudinal history
data and a fitted `mjoint` object.

## Usage

``` r
dynSurv(
  object,
  newdata,
  newSurvData = NULL,
  u = NULL,
  horizon = NULL,
  type = "first-order",
  M = 200,
  scale = 2,
  ci,
  progress = TRUE
)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- newdata:

  a list of `data.frame` objects for each longitudinal outcome for a
  single new patient in which to interpret the variables named in the
  `formLongFixed` and `formLongRandom` formulae of `object`. As per
  [`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
  the `list` structure enables one to include multiple longitudinal
  outcomes with different measurement protocols. If the multiple
  longitudinal outcomes are measured at the same time points for each
  patient, then a `data.frame` object can be given instead of a `list`.
  It is assumed that each data frame is in long format.

- newSurvData:

  a `data.frame` in which to interpret the variables named in the
  `formSurv` formulae from the `mjoint` object. This is optional, and if
  omitted, the data will be searched for in `newdata`. Note that no
  event time or censoring indicator data are required for dynamic
  prediction. Defaults to `newSurvData=NULL`.

- u:

  an optional time that must be greater than the last observed
  measurement time. If omitted (default is `u=NULL`), then conditional
  failure probabilities are reported for *all* observed failure times in
  the `mjoint` object data from the last known follow-up time of the
  subject.

- horizon:

  an optional horizon time. Instead of specifying a specific time `u`
  relative to the time origin, one can specify a horizon time that is
  relative to the last known follow-up time. The prediction time is
  essentially equivalent to `horizon` + \\t\_{obs}\\, where \\t\_{obs}\\
  is the last known follow-up time where the patient had not yet
  experienced the event. Default is `horizon=NULL`. If `horizon` is
  non-`NULL`, then the output will be reported still in terms of
  absolute time (from origin), `u`.

- type:

  a character string for whether a first-order (`type="first-order"`) or
  Monte Carlo simulation approach (`type="simulated"`) should be used
  for the dynamic prediction. Defaults to the computationally faster
  first-order prediction method.

- M:

  for `type="simulated"`, the number of simulations to performs. Default
  is `M=200`.

- scale:

  a numeric scalar that scales the variance parameter of the proposal
  distribution for the Metropolis-Hastings algorithm, which therefore
  controls the acceptance rate of the sampling algorithm.

- ci:

  a numeric value with value in the interval \\(0, 1)\\ specifying the
  confidence interval level for predictions of `type='simulated'`. If
  missing, defaults to `ci=0.95` for a 95% confidence interval. If
  `type='first-order'` is used, then this argument is ignored.

- progress:

  logical: should a progress bar be shown on the console to indicate the
  percentage of simulations completed? Default is `progress=TRUE`.

## Value

A list object inheriting from class `dynSurv`. The list returns the
arguments of the function and a `data.frame` with first column (named
`u`) denoting times and the subsequent columns returning summary
statistics for the conditional failure probabilities For
`type="first-order"`, a single column named `surv` is appended. For
`type="simulated"`, four columns named `mean`, `median`, `lower` and
`upper` are appended, denoting the mean, median and lower and upper
confidence intervals from the Monte Carlo draws. Additional objects are
returned that are used in the intermediate calculations.

## Details

Dynamic predictions for the time-to-event data sub-model based on an
observed measurement history for the longitudinal outcomes of a new
subject are based on either a first-order approximation or Monte Carlo
simulation approach, both of which are described in Rizopoulos (2011).
Namely, given that the subject was last observed at time *t*, we
calculate the conditional survival probability at time \\u \> t\\ as

\$\$P\[T \ge u \| T \ge t; y, \theta\] \approx \frac{S(u \| \hat{b};
\theta)}{S(t \| \hat{b}; \theta)},\$\$

where \\T\\ is the failure time for the new subject, \\y\\ is the
stacked-vector of longitudinal measurements up to time *t* and \\S(u \|
\hat{b}; \theta)\\ is the survival function.

**First order predictions**

For `type="first-order"`, \\\hat{b}\\ is the mode of the posterior
distribution of the random effects given by

\$\$\hat{b} = {\arg \max}\_b f(b \| y, T \ge t; \theta).\$\$

The predictions are based on plugging in \\\theta = \hat{\theta}\\,
which is extracted from the `mjoint` object.

**Monte Carlo simulation predictions**

For `type="simulated"`, \\\theta\\ is drawn from a multivariate normal
distribution with means \\\hat{\theta}\\ and variance-covariance matrix
both extracted from the fitted `mjoint` object via the
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) functions. \\\hat{b}\\ is
drawn from the the posterior distribution of the random effects

\$\$f(b \| y, T \ge t; \theta)\$\$

by means of a Metropolis-Hasting algorithm with independent multivariate
non-central *t*-distribution proposal distributions with non-centrality
parameter \\\hat{b}\\ from the first-order prediction and
variance-covariance matrix equal to `scale` \\\times\\ the inverse of
the negative Hessian of the posterior distribution. The choice of
`scale` can be used to tune the acceptance rate of the
Metropolis-Hastings sampler. This simulation algorithm is iterated `M`
times, at each time calculating the conditional survival probability.

## References

Rizopoulos D. Dynamic predictions and prospective accuracy in joint
models for longitudinal and time-to-event data. *Biometrics*. 2011;
**67**: 819–829.

Taylor JMG, Park Y, Ankerst DP, Proust-Lima C, Williams S, Kestin L, et
al. Real-time individual predictions of prostate cancer recurrence using
joint models. *Biometrics*. 2013; **69**: 206–13.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`dynLong`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md),
and
[`plot.dynSurv`](https://graemeleehickey.github.io/joineRML/reference/plot.dynSurv.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a joint model with bivariate longitudinal outcomes

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

fit2 <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    verbose = TRUE)

hvd2 <- droplevels(hvd[hvd$num == 1, ])
dynSurv(fit2, hvd2)
dynSurv(fit2, hvd2, u = 7) # survival at 7-years only

out <- dynSurv(fit2, hvd2, type = "simulated")
out
} # }
```
