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
  progress = TRUE,
  max.proposals = 10000
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
  survival probabilities are reported for *all* observed failure times
  in the `mjoint` object data from the last known follow-up time of the
  subject.

- horizon:

  an optional horizon time. Instead of specifying a specific time `u`
  relative to the time origin, one can specify a horizon time that is
  relative to the last known follow-up time. The prediction time is
  essentially equivalent to `horizon` + \\t\_{obs}\\, where \\t\_{obs}\\
  is the last known follow-up time where the patient had not yet
  experienced the event. Default is `horizon=NULL`. If `horizon` is
  non-`NULL`, then the output will be reported still in terms of
  absolute time (from origin), `u`. Supplied horizons must be numeric,
  finite, and strictly positive.

- type:

  a character string for whether a first-order (`type="first-order"`) or
  Monte Carlo simulation approach (`type="simulated"`) should be used
  for the dynamic prediction. Defaults to the computationally faster
  first-order prediction method.

- M:

  for `type="simulated"`, the number of simulations. Must be a positive
  integer. Default is `M=200`.

- scale:

  deprecated; ignored. Supplying it produces a warning. The conditional
  rejection sampler does not require a proposal scale.

- ci:

  a numeric value with value in the interval \\(0, 1)\\ specifying the
  confidence interval level for predictions of `type='simulated'`. If
  missing, defaults to `ci=0.95` for a 95% confidence interval. If
  `type='first-order'` is used, then this argument is ignored.

- progress:

  logical: should a progress bar be shown on the console to indicate the
  percentage of simulations completed? Default is `progress=TRUE`.

- max.proposals:

  maximum number of proposals per parameter draw for `type="simulated"`.
  Must be a positive integer; defaults to `10000`. The prediction stops
  if this limit is exhausted.

## Value

A list object inheriting from class `dynSurv`. The list returns the
arguments of the function and a `data.frame` with first column (named
`u`) denoting times and the subsequent columns returning summary
statistics for the conditional survival probabilities. For
`type="first-order"`, a single column named `surv` is appended. For
`type="simulated"`, four columns named `mean`, `median`, `lower` and
`upper` are appended, denoting the mean, median and lower and upper
confidence limits from the Monte Carlo draws. Additional objects are
returned that are used in the intermediate calculations. For simulated
predictions, `proposals` records the number of proposals for each of the
`M` accepted random-effects draws, `accept` is `M / sum(proposals)`, and
`sampler` is `"conditional-rejection"`. The former `scale` setting is
ignored. The `b.hat` component contains the conditional mode for
first-order predictions and is `NULL` for simulated predictions.

## Details

Dynamic predictions for the time-to-event data sub-model based on an
observed measurement history for the longitudinal outcomes of a new
subject are based on either a first-order approximation or Monte Carlo
simulation approach. See Rizopoulos (2011) for the prediction framework.
Given that the subject was last observed at time *t*, we calculate the
conditional survival probability at time \\u \> t\\ as

\$\$P\[T \ge u \| T \ge t; y, \theta\] \approx \frac{S(u \| \hat{b};
\theta)}{S(t \| \hat{b}; \theta)},\$\$

where \\T\\ is the failure time for the new subject, \\y\\ is the
stacked-vector of longitudinal measurements up to time *t* and \\S(u \|
\hat{b}; \theta)\\ is the survival function. The ratio is evaluated from
hazard increments between \\t\\ and \\u\\, avoiding division of survival
probabilities that may both round to zero.

Prediction matrices use the transformations, factor levels, and
contrasts retained when the model was fitted. Changing
`options("contrasts")` does not change predictions from that fit. New
factor levels are rejected. Older fits are checked against their stored
random-effects matrices; inconsistent bases at measurement and event
times require refitting with
[`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md).
Transformations that depend on the data must support R's
[`makepredictcall`](https://rdrr.io/r/stats/makepredictcall.html)
mechanism or use fixed constants.

**First order predictions**

For `type="first-order"`, \\\hat{b}\\ is the mode of the posterior
distribution of the random effects given by

\$\$\hat{b} = {\arg \max}\_b f(b \| y, T \ge t; \theta).\$\$

The predictions are based on plugging in \\\theta = \hat{\theta}\\,
which is extracted from the `mjoint` object.

**Monte Carlo simulation predictions**

For `type="simulated"`, the non-hazard parameters are drawn jointly from
a multivariate normal approximation with residual variances on the log
scale. Write \\v_k = \sigma_k^2\\ and let \\g(\theta)\\ replace each
\\v_k\\ by \\\log(v_k)\\, leaving the other parameters unchanged. The
mean is \\g(\hat{\theta})\\, and the covariance matrix is \\J V J^T\\,
where \\V\\ is returned by `vcov(object)` and \\J\\ is diagonal, with
entries \\1 / \hat{v}\_k\\ for the residual variances and 1 elsewhere.
This delta-method transformation includes the covariances with the other
parameters. Exponentiating the log-variance draws gives positive
residual variances, with medians equal to their fitted estimates. Their
original-scale means and covariances are not preserved exactly. Fitted
estimates and reported standard errors retain their original scales.
Invalid fitted variances, unusable covariance matrices, and numerical
overflow or underflow in the draws produce an error.

This is an asymptotic approximation; positivity alone does not establish
its accuracy in small samples or near a variance boundary.
Random-effects covariance matrices and baseline-hazard increments retain
their existing adjustments to ensure positive definiteness and
positivity, respectively.

For each parameter draw, random effects are sampled from \\f(b \| y, T
\ge t; \theta)\\ by Gaussian rejection sampling. The longitudinal
conditional distribution is \\N(m,V)\\. Its density is multiplied by the
survival contribution \\\exp\\\ell(b)\\\\, where \\\ell(b) = -\sum_j h_j
\exp(a_j^T b)\\ is concave and \\h_j \ge 0\\. At an approximate
conditional mode \\b_0\\, put \\g_0 = \nabla\ell(b_0)\\. A tangent bound
gives the proposal \\N(m + V g_0,V)\\, accepted with log probability
\\\ell(b) - \ell(b_0) - g_0^T(b-b_0)\\. This samples the conditional
distribution without burn-in; an approximate mode affects efficiency,
rather than the target distribution. Acceptance is evaluated on the log
scale to avoid cancellation and overflow.

Each parameter draw is held fixed until a random-effects draw is
accepted. Exhausting `max.proposals` stops the prediction and identifies
the simulation and conditioning time; parameter draws are never
discarded and replaced. Strong survival conditioning, particularly in
several independent random-effects directions, can require many
proposals. Increasing `max.proposals` allows more attempts but does not
improve the approximation used to sample the model parameters. The
conditional sampler alone does not establish prediction-interval
coverage.

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
