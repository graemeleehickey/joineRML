# joineRML

## Introduction

The `joineRML` package implements methods for analyzing data from
*multiple* longitudinal studies in which the responses from each subject
consists of time-sequences of repeated measurements and a possibly
censored time-to-event outcome. The modelling framework for the repeated
measurements is the multivariate linear mixed effects model. The model
for the time-to-event outcome is a Cox proportional hazards model with
log-Gaussian frailty. Stochastic dependence is captured by allowing the
Gaussian random effects of the linear model to be correlated with the
frailty term of the Cox proportional hazards model. For full details of
the model, please consult the technical vignette by running

``` r

vignette("technical", package = "joineRML")
```

## Heart valve data

### Data

The simplest way to explain the concepts of the package is through an
example. `joineRML` comes with the data set `heart.valve`. Details of
this data can be found in the help file by running the command

``` r

help("heart.valve", package = "joineRML")
```

This data is in so-called *long* or *unbalanced* format:

``` r

library("joineRML")
data("heart.valve")
head(heart.valve)
#>   num sex      age      time    fuyrs status grad log.grad   lvmi log.lvmi ef
#> 1   1   0 75.06027 0.0109589 4.956164      0   10 2.302585 118.98 4.778955 93
#> 2   1   0 75.06027 3.6794520 4.956164      0   10 2.302585 118.98 4.778955 93
#> 3   1   0 75.06027 4.6958900 4.956164      0   10 2.302585 137.63 4.924569 93
#> 4   2   0 45.79452 6.3643840 9.663014      0   14 2.639057 114.93 4.744323 68
#> 5   2   0 45.79452 7.3041100 9.663014      0    9 2.197225 109.80 4.698661 70
#> 6   2   0 45.79452 8.3013700 9.663014      0   12 2.484907 157.40 5.058790 56
#>    bsa lvh prenyha redo size con.cabg creat dm acei lv emergenc hc sten.reg.mix
#> 1 1.77   1       3    0   27        1   103  0    1  1        0  0            1
#> 2 1.77   1       3    0   27        1   103  0    1  1        0  0            1
#> 3 1.77   1       3    0   27        1   103  0    1  1        0  0            1
#> 4 1.92   1       1    1   22        0    76  0    0  2        0  0            1
#> 5 1.92   1       1    1   22        0    76  0    0  2        0  0            1
#> 6 1.92   1       1    1   22        0    76  0    0  2        0  0            1
#>                hs
#> 1 Stentless valve
#> 2 Stentless valve
#> 3 Stentless valve
#> 4       Homograft
#> 5       Homograft
#> 6       Homograft
```

The data refer to 256 patients and are stored in the unbalanced format,
which is convenient here because measurement times were unique to each
subject. The data are stored as a single R object, `heart.valve`, which
is a data frame of dimension 988 by 25. The average number of repeated
measurements per subject is therefore 988/256 = 3.86. As with any
unbalanced data set, values of time-constant variables are repeated over
all rows that refer to the same subject. The dimensionality of the data
set can be confirmed by a call to the
[`dim()`](https://rdrr.io/r/base/dim.html) function, whilst the names of
the 25 variables can be listed by a call to the
[`names()`](https://rdrr.io/r/base/names.html) function:

``` r

dim(heart.valve)
#> [1] 988  25
names(heart.valve)
#>  [1] "num"          "sex"          "age"          "time"         "fuyrs"       
#>  [6] "status"       "grad"         "log.grad"     "lvmi"         "log.lvmi"    
#> [11] "ef"           "bsa"          "lvh"          "prenyha"      "redo"        
#> [16] "size"         "con.cabg"     "creat"        "dm"           "acei"        
#> [21] "lv"           "emergenc"     "hc"           "sten.reg.mix" "hs"
```

We will only analyse a subset of this data, namely records with
case-complete data for heart valve gradient (`grad`) and left
ventricular mass index (`lvmi`):

``` r

hvd <- heart.valve[!is.na(heart.valve$grad) & !is.na(heart.valve$lvmi), ]
```

Strictly speaking, this is not necessary because `joineRML` can handle
the situation of different measurement schedules *within* subjects That
is, a subject does not need to have all multiple longitudinal outcomes
recorded at each visit. It is conceivable that some biomarkers will be
measured more or less frequently than others. For example, invasive
measurements may only be recorded annually, whereas a simple biomarker
measurement might be recorded more frequently. `joineRML` can handle
this situation by specifying each longitudinal outcome its own data
frame.

Further to that, we only select the first 50 individuals to speed up
these examples:

``` r

hvd <- hvd[hvd$num <= 50, ]
```

### Model fitting

The main function in the `joineRML` package is the
[`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
function. Its main (required) arguments are:

- `formLongFixed`: a list (of length equal to the number of longitudinal
  outcome types considered) of two-sided formulae specifying the
  response on the left-hand side and the mean linear predictor terms for
  the fixed effects in the linear mixed models on the right-hand side.

- `formLongRandom`: a list (of same length as `formLongFixed`) of
  one-sided formulae specifying the model for random effects in the
  linear mixed models.

- `formSurv`: a formula specifying baseline covariates in the
  proportional hazards regression model for the time-to-event outcome,
  with one common baseline hazard. The latent association structure is
  added separately. Factors, interactions, and unpenalised
  transformations are supported. Cox-specific terms such as
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html),
  [`offset()`](https://rdrr.io/r/stats/offset.html),
  [`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html), `tt()`,
  [`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html),
  [`pspline()`](https://rdrr.io/pkg/survival/man/pspline.html), and
  [`ridge()`](https://rdrr.io/pkg/survival/man/ridge.html) are not
  supported and are rejected before fitting, including calls prefixed
  with `survival::` or `stats::`.

- `data`: a list (of same length as `formLongFixed`) of data.frames; one
  for each longitudinal outcome. It is assumed that the event time data
  is in the first data.frame (i.e. `data[[1]]`), unless the argument
  `survData` (which defaults to `NULL`) is specified. If $`K>1`$ and the
  data are balanced within patients (i.e. multiple markers measured at
  common measurement times), then one can specify `data` as a data frame
  rather than as a list.

- `timeVar`: the column name indicating the time variable in the linear
  mixed effects model. If $`K>1`$ and the data frames have different
  column names for time, then `timeVar` can alternatively be specified
  as a vector of strings of length $`K`$.

Subject IDs must be non-missing and identify the same subjects in every
dataset. Records are matched by their ID labels, so row order and
factor-level order can differ between datasets. Measurement schedules
may differ between longitudinal outcomes, but each subject must
contribute observations for every outcome. A separate `survData` must
contain exactly one record per subject, with the same IDs as the
longitudinal data and no missing values in the variables used by the
survival model. At least one observed event is required.

We can fit a bivariate joint model to the log-transformed valve gradient
and LVMI indices in the `hvd` subset using

``` r

fit <- mjoint(
  formLongFixed = list("grad" = log.grad ~ time + sex + hs, 
                       "lvmi" = log.lvmi ~ time + sex),
  formLongRandom = list("grad" = ~ 1 | num,
                        "lvmi" = ~ time | num),
  formSurv = Surv(fuyrs, status) ~ age,
  data = list(hvd, hvd),
  inits = list("gamma" = c(0.11, 1.51, 0.80)),
  timeVar = "time",
  seed = 12345)
#> Running multivariate LMM EM algorithm to establish initial parameters...
#> Finished multivariate LMM EM algorithm...
#> EM algorithm has converged!
#> Calculating post model fit statistics...
```

The `seed` argument makes the whole fit reproducible, including post-fit
calculations, for the same data, settings, random-number generator, and
software versions. It restores the previous R random-number state when
the fit finishes. The default `seed = NULL` uses the current state as
usual. The supplied seed does not control subsequent calls to
[`bootSE()`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md),
[`dynLong()`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md),
or
[`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md);
set the seed separately before those calculations. See
[`help("mjoint")`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
for details of random-number generation.

Details on the model estimation algorithm are provided in the technical
details vignette. We note here that this is not necessarily the most
appropriate model for the data, and is included only for the purposes of
demonstration. There are a number of other useful arguments in the
`mjoint` function; for example, `inits` for specifying (partial) initial
values, `control` for controlling the optimization algorithm, and
`verbose` for monitoring the convergence output in real-time. A full
list of all arguments with explanation are given in the help
documentation, accessed by running
[`help("mjoint")`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md).

### Post-fit analysis

Once we have a fitted `mjoint` object, we can begin to extract relevant
information from it. Most summary statistics are available from the
`summary` function:

``` r

summary(fit)
#> 
#> Call:
#> mjoint(formLongFixed = list(grad = log.grad ~ time + sex + hs, 
#>     lvmi = log.lvmi ~ time + sex), formLongRandom = list(grad = ~1 | 
#>     num, lvmi = ~time | num), formSurv = Surv(fuyrs, status) ~ 
#>     age, data = list(hvd, hvd), timeVar = "time", inits = list(gamma = c(0.11, 
#>     1.51, 0.8)), seed = 12345L)
#> 
#> Data Descriptives:
#> 
#> Event Process
#>     Number of subjects: 43 
#>     Number of events: 8 (18.6%)
#> 
#> Longitudinal Process
#>     Number of longitudinal outcomes: K = 2 
#>     Number of observations:
#>       Outcome 1 (grad): n = 144
#>       Outcome 2 (lvmi): n = 144
#> 
#> Joint Model Summary:
#> 
#> Longitudinal Process: Multivariate linear mixed-effects model
#>      log.grad ~ time + sex + hs, random = ~1 | num
#>      log.lvmi ~ time + sex, random = ~time | num
#> Event Process: Cox proportional hazards model
#>      Surv(fuyrs, status) ~ age
#> Model fit statistics:
#>    log.Lik      AIC      BIC
#>  -153.2218 342.4435 374.1451
#> 
#> Variance Components:
#> 
#> Random effects variance covariance matrix
#>               (Intercept)_1 (Intercept)_2     time_2
#> (Intercept)_1     0.0958880     0.0596100 -0.0051222
#> (Intercept)_2     0.0596100     0.1471500 -0.0049885
#> time_2           -0.0051222    -0.0049885  0.0011174
#>   Standard Deviations: 0.30966 0.3836 0.033427 
#> 
#> Residual standard deviations:
#>      grad      lvmi 
#> 0.5171520 0.1581865 
#> 
#> Coefficient Estimates:
#> 
#> Longitudinal sub-model:
#>                       Value Std.Err z-value p-value
#> (Intercept)_1        2.5219  0.1654 15.2481 <0.0001
#> time_1               0.0095  0.0307  0.3102  0.7564
#> sex_1                0.0344  0.2312  0.1489  0.8817
#> hsStentless valve_1  0.0825  0.1990  0.4146  0.6784
#> (Intercept)_2        4.9929  0.0818 61.0612 <0.0001
#> time_2               0.0258  0.0129  1.9974  0.0458
#> sex_2               -0.1967  0.2105 -0.9344  0.3501
#> 
#> Time-to-event sub-model:
#>          Value Std.Err z-value p-value
#> age     0.1366  0.0665  2.0527  0.0401
#> gamma_1 2.8472  3.7330  0.7627  0.4456
#> gamma_2 3.6801  3.4488  1.0671  0.2859
#> 
#> Algorithm Summary:
#>     Total computational time: 4.1 secs 
#>     EM algorithm computational time: 3.6 secs 
#>     Convergence status: converged
#>     Convergence criterion: sas 
#>     Final Monte Carlo sample size: 3674 
#>     Standard errors calculated using method: approx
```

The residual quantities printed by `summary(fit)` are standard
deviations, also returned by `sigma(fit)`; the stored `coef(fit)$sigma2`
values are variances.

One can also extract the coefficients, fixed effects, and random effects
using standard generic functions:

``` r

coef(fit)
#> $D
#>               (Intercept)_1 (Intercept)_2       time_2
#> (Intercept)_1    0.09588800   0.059610212 -0.005122220
#> (Intercept)_2    0.05961021   0.147152034 -0.004988491
#> time_2          -0.00512222  -0.004988491  0.001117369
#> 
#> $beta
#>       (Intercept)_1              time_1               sex_1 hsStentless valve_1 
#>         2.521878812         0.009511296         0.034415375         0.082506057 
#>       (Intercept)_2              time_2               sex_2 
#>         4.992936048         0.025845910        -0.196713101 
#> 
#> $sigma2
#>   sigma2_1   sigma2_2 
#> 0.26744623 0.02502296 
#> 
#> $haz
#> [1] 0.002239315 0.002463859 0.005677867 0.006542868 0.008955061 0.011749720
#> [7] 0.020522447 0.266416490
#> 
#> $gamma
#>       age   gamma_1   gamma_2 
#> 0.1365893 2.8471790 3.6800986
fixef(fit, process = "Longitudinal")
#>       (Intercept)_1              time_1               sex_1 hsStentless valve_1 
#>         2.521878812         0.009511296         0.034415375         0.082506057 
#>       (Intercept)_2              time_2               sex_2 
#>         4.992936048         0.025845910        -0.196713101
fixef(fit, process = "Event")
#>       age   gamma_1   gamma_2 
#> 0.1365893 2.8471790 3.6800986
head(ranef(fit))
#>    (Intercept)_1 (Intercept)_2       time_2
#> 1    -0.18893011   -0.24357531  0.005770637
#> 10   -0.16257082   -0.04703603  0.003013689
#> 11    0.04615938    0.05514878  0.009264502
#> 12    0.49281848    0.36760602 -0.028566714
#> 13    0.36964501    0.32790364 -0.035443210
#> 14   -0.15438221   -0.56119773  0.018966600
```

Although a model fit may indicate convergence, it is generally a good
idea to examine the convergence plots. These can be viewed using the
`plot` function for each group of model parameters.

``` r

plot(fit, params = "gamma")
```

![](joineRML_files/figure-html/hvd_model_conv-1.png)

``` r

plot(fit, params = "beta")
```

![](joineRML_files/figure-html/hvd_model_conv-2.png)

### Dynamic predictions

Use
[`dynLong()`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md)
for future longitudinal values and
[`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
for conditional survival probabilities, given a subject’s observed
measurement history. Their default `type = "first-order"` uses the
fitted parameters and the conditional random-effects mode. Setting
`type = "simulated"` also samples parameters and draws random effects
conditional on each sampled parameter vector; set the seed separately
before each prediction call for reproducibility. See
[`help("dynLong")`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md)
and
[`help("dynSurv")`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
for examples and arguments.

Prediction uses the transformations and factor coding established during
fitting. Random-effects time terms such as `poly(time, 2)`,
`splines::ns(time, df = 2)`, and `scale(time)` retain their fitted
coefficients, knots, or scaling constants at event and prediction times,
including when predicting at a single time. Changing
`options("contrasts")` after fitting does not change predictions, and
new factor levels are rejected. Custom transformations that depend on
the data must support R’s
[`makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
mechanism or use fixed constants; an expression such as
`time - mean(time)` recalculates its centring whenever it is evaluated.
Older fits are checked against their stored random-effects matrices. If
measurement-time and event-time bases disagree, the model must be
refitted with
[`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md);
recalculating predictions alone cannot correct the fitted model.

For simulated predictions, residual variances are sampled jointly with
the other non-hazard parameters on the log-variance scale. The full
covariance matrix is transformed using the delta method, including
covariances between residual variances and other parameters, and the
log-variance draws are then exponentiated. This gives positive residual
variances with medians equal to the fitted estimates; it does not
preserve their original-scale means or covariances exactly. Fitted
estimates, [`vcov()`](https://rdrr.io/r/stats/vcov.html), and reported
standard errors retain their existing scales. Invalid fitted variances,
unusable covariance matrices, and numerical overflow or underflow
produce an error. The approximation remains asymptotic, so positive
variance draws alone do not establish reliable prediction intervals in
small samples or near a variance boundary.

For each sampled parameter vector, Gaussian rejection sampling draws
random effects conditional on the observed measurements and survival to
the last observation time. The parameter vector stays fixed until a
random-effects draw is accepted, so no burn-in is needed. The proposal
uses a tangent bound to the concave log-survival contribution; an
approximate conditional mode affects efficiency rather than the
conditional distribution being sampled. The former `scale` argument is
deprecated and ignored, with a warning if supplied.

The default `max.proposals = 10000` limits attempts for each parameter
draw. If the limit is exhausted, prediction stops and reports the
simulation number and conditioning time; it does not discard that
parameter draw. Strong survival conditioning can require many proposals.
Check the fit and its uncertainty before increasing the limit. The
returned `proposals` vector gives the number of attempts for each
accepted draw, and `accept` is `M / sum(proposals)`. These describe
sampling efficiency, not the accuracy of the parameter approximation or
prediction-interval coverage. Conditional survival probabilities are
calculated from hazard increments between the conditioning and
prediction times, avoiding division of probabilities that may both round
to zero.

### Comparing fitted models

Use [`anova()`](https://rdrr.io/r/stats/anova.html) to compare nested
covariate models fitted to the same subjects, longitudinal responses and
measurement times, and survival responses. The response transformations,
random-effects structure, and association structure must also agree.
Here, the reduced model omits `hs` from the valve-gradient fixed effects
and retains all other terms. The existing `fit` is the full model.
Setting `pfs = FALSE` skips the reduced model’s post-fit standard errors
and random-effect predictions, which are not needed for this comparison.

``` r

fit.reduced <- update(fit,
  formLongFixed = list("grad" = log.grad ~ time + sex,
                       "lvmi" = log.lvmi ~ time + sex),
  pfs = FALSE,
  seed = 23456)
#> Running multivariate LMM EM algorithm to establish initial parameters...
#> Finished multivariate LMM EM algorithm...
#> EM algorithm has converged!
```

Supply models from reduced to full. Both must have converged;
[`anova()`](https://rdrr.io/r/stats/anova.html) also checks that the
observations and model structures permit the comparison.

``` r

comparison <- anova(fit.reduced, fit)
knitr::kable(comparison, digits = 4)
```

|             |  df |    logLik |      AIC |      BIC | L.Ratio |  Df | Pr(\>Chisq) | LR.change |
|:------------|----:|----------:|---------:|---------:|--------:|----:|------------:|----------:|
| fit.reduced |  17 | -153.4123 | 340.8246 | 370.7650 |      NA |  NA |          NA |        NA |
| fit         |  18 | -153.2311 | 342.4622 | 374.1638 |  0.3624 |   1 |      0.5472 |     2e-04 |

The test uses an approximate chi-squared reference distribution, with
degrees of freedom equal to the number of additional identifiable
coefficients. The joint likelihoods are re-evaluated at the fitted
coefficients using common deterministic Halton integration points,
without refitting either model. Consequently, the table’s
log-likelihoods, AIC, and BIC may differ from the values stored in the
fitted objects. BIC uses the number of subjects, not the number of
longitudinal measurements.

`LR.change` is the absolute change in the likelihood-ratio statistic
when the integration size increases from half of `nMC` to `nMC` (default
10,000). A negative statistic or a change exceeding `tol` (default 0.1)
produces a warning and an `NA` p-value. To investigate integration
error, increase `nMC`, for example:

``` r

anova(fit.reduced, fit, nMC = 20000)
```

The tolerance is a numerical diagnostic, not a significance threshold or
a guaranteed error bound. A negative statistic can also indicate that
the fits need tighter convergence settings or better initial values.
Passing the diagnostic does not establish that the maximum was reached
or that the asymptotic approximation is adequate. This small dataset is
used to demonstrate the workflow, rather than to draw clinical
conclusions.

For a descriptive comparison using the stored likelihoods, without
further integration or a likelihood-ratio test, use:

``` r

anova(fit.reduced, fit, test = FALSE)
#> Joint model comparison (maximum likelihood)
#> 
#> Stored likelihoods; no likelihood-ratio tests.
#> 
#>             df  logLik    AIC    BIC
#> fit.reduced 17 -153.43 340.87 370.81
#> fit         18 -153.22 342.44 374.15
```

This option also permits non-nested models or different random-effects
structures, provided the models have converged and use the same
observations and response transformations. Likelihood-ratio tests of
changes to the random-effects structure, including zero variance
components, are not supported. The [tidying
vignette](https://graemeleehickey.github.io/joineRML/articles/joineRML-tidy.md)
shows how to extract model components and format coefficient and
model-comparison tables.

### Bootstrap standard errors

Once an `mjoint` model has converged, and assuming the `pfs` argument is
`TRUE` (default), then approximated standard errors are calculated based
on the empirical information matrix of the profile likelihood at the
maximizer. Theoretically, these standard errors will be underestimated
(see the technical vignette). In principle, residual Monte Carlo error
will oppose this through an increase in uncertainty.

``` r

fit.se <- bootSE(fit,
                 nboot = 100,
                 ncores = 1L)
```

Bootstrapping is a computationally intensive method, possibly taking
many hours to fit. For this reason, one can relax the control parameter
constraints on the optimization algorithm for each bootstrap model;
however, this will be at the possible expense of inflated standard
errors due to Monte Carlo error.

We can call the `bootSE` object to interrogate it

``` r

fit.se
```

The residual estimates and percentile confidence intervals printed by
`fit.se` are on the variance scale; their `SE` column gives standard
errors of the variance estimates. These are distinct from the residual
standard deviations printed by `summary(fit)`.

Alternatively, re-run the `summary` command, passing the additional
argument of `bootSE = fit.se`:

``` r

summary(fit, bootSE = fit.se)
```

## Univariate joint models: `joineRML` versus `joineR`

There are a growing number of software options for fitting joint models
of a single longitudinal outcome and a single time-to-event outcome;
what we call here *univariate* joint models. `joineR` (version 1.2.7) is
one package available in R for fitting such models, however `joineRML`
can fit these models too, since the univariate model is simply a special
case of the multivariate model. It is useful to contrast these two
packages. There are theoretical and practical implementation differences
between the packages beyond just univariate versus multivariate
capability:

- `joineR` uses Gauss-Hermite quadrature (with 3 nodes) for numerical
  integration, whereas `joineRML` uses Monte Carlo integration with an
  automated selection of sample size.

- `joineR` only allows for random-intercept models, random-intercept and
  random-slope models, or a quadratic model. `joineRML`, on the other
  hand, allows for any random effects structure.

- `joineR` only allows for specification of convergence based on an
  absolute difference criterion.

- `joineR` does not calculate approximate standard errors, and instead
  requires a bootstrap approach be used after the model fit.

- The current version of `joineR` requires a data pre-processing step in
  order to generate a `joint.data` object, whereas `joineRML` can work
  straight from the data frame.

To fit a univariate model in `joineR` we run the following code for the
`hvd` data

``` r

library(joineR, quietly = TRUE)

hvd.surv <- UniqueVariables(hvd, var.col = c("fuyrs", "status"), id.col = "num")
hvd.cov <- UniqueVariables(hvd, "age", id.col = "num")
hvd.long <- hvd[, c("num", "time", "log.lvmi")]

hvd.jd <- jointdata(longitudinal = hvd.long, 
                    baseline = hvd.cov, 
                    survival = hvd.surv, 
                    id.col = "num", 
                    time.col = "time")

fit.joiner <- joint(data = hvd.jd,
                    long.formula = log.lvmi ~ time + age, 
                    surv.formula = Surv(fuyrs, status) ~ age, 
                    model = "intslope")

summary(fit.joiner)
#> 
#> Call:
#> joint(data = hvd.jd, long.formula = log.lvmi ~ time + age, surv.formula = Surv(fuyrs, 
#>     status) ~ age, model = "intslope")
#> 
#> Random effects joint model
#>  Data: hvd.jd 
#>  Log-likelihood: -30.44117 
#> 
#> Longitudinal sub-model fixed effects: log.lvmi ~ time + age                       
#> (Intercept) 4.601687012
#> time        0.027065527
#> age         0.005898348
#> 
#> Survival sub-model fixed effects: Surv(fuyrs, status) ~ age             
#> age 0.1192532
#> 
#> Latent association:                
#> gamma_0 3.496122
#> 
#> Variance components:
#>          U_0          U_1     Residual 
#> 0.1358091795 0.0008925357 0.0254863247 
#> 
#> Convergence at iteration: 43 
#> 
#> Number of observations: 144 
#> Number of groups: 43
```

To fit a univariate model in `joineRML` we run the following code for
the `hvd` data

``` r

fit.joinerml <- mjoint(formLongFixed = log.lvmi ~ time + age,
                       formLongRandom = ~ time | num,
                       formSurv = Surv(fuyrs, status) ~ age,
                       data = hvd,
                       timeVar = "time",
                       seed = 123)
#> EM algorithm has converged!
#> Calculating post model fit statistics...

summary(fit.joinerml)
#> 
#> Call:
#> mjoint(formLongFixed = log.lvmi ~ time + age, formLongRandom = ~time | 
#>     num, formSurv = Surv(fuyrs, status) ~ age, data = hvd, timeVar = "time", 
#>     seed = 123L)
#> 
#> Data Descriptives:
#> 
#> Event Process
#>     Number of subjects: 43 
#>     Number of events: 8 (18.6%)
#> 
#> Longitudinal Process
#>     Number of longitudinal outcomes: K = 1 
#>     Number of observations:
#>       Outcome 1: n = 144
#> 
#> Joint Model Summary:
#> 
#> Longitudinal Process: Univariate linear mixed-effects model
#>      log.lvmi ~ time + age, random = ~time | num
#> Event Process: Cox proportional hazards model
#>      Surv(fuyrs, status) ~ age
#> Model fit statistics:
#>    log.Lik      AIC      BIC
#>  -30.42324 78.84648 94.69729
#> 
#> Variance Components:
#> 
#> Random effects variance covariance matrix
#>               (Intercept)_1      time_1
#> (Intercept)_1     0.1357700 -0.00364490
#> time_1           -0.0036449  0.00087096
#>   Standard Deviations: 0.36847 0.029512 
#> 
#> Residual standard deviations:
#>   sigma_1 
#> 0.1598896 
#> 
#> Coefficient Estimates:
#> 
#> Longitudinal sub-model:
#>                Value Std.Err z-value p-value
#> (Intercept)_1 4.6021  0.4178 11.0163 <0.0001
#> time_1        0.0271  0.0105  2.5743  0.0100
#> age_1         0.0059  0.0066  0.8982  0.3691
#> 
#> Time-to-event sub-model:
#>          Value Std.Err z-value p-value
#> age     0.1190  0.0462  2.5761  0.0100
#> gamma_1 3.4823  1.6559  2.1030  0.0355
#> 
#> Algorithm Summary:
#>     Total computational time: 1.7 secs 
#>     EM algorithm computational time: 1.5 secs 
#>     Convergence status: converged
#>     Convergence criterion: sas 
#>     Final Monte Carlo sample size: 3784 
#>     Standard errors calculated using method: approx
```

In addition to just comparing model parameter estimates, we can also
extract the predicted (or posterior) random effects from each model and
plot them.

``` r

id <- as.numeric(row.names(fit.joiner$coefficients$random))
id.ord <- order(id) # joineR rearranges patient ordering during EM fit
par(mfrow = c(1, 2))
plot(fit.joiner$coefficients$random[id.ord, 1], ranef(fit.joinerml)[, 1],
     main = "Predicted random intercepts",
     xlab = "joineR", ylab = "joineRML")
grid()
abline(a = 0, b = 1, col = 2, lty = "dashed")
plot(fit.joiner$coefficients$random[id.ord, 2], ranef(fit.joinerml)[, 2],
     main = "Predicted random slopes",
     xlab = "joineR", ylab = "joineRML")
grid()
abline(a = 0, b = 1, col = 2, lty = "dashed")
```

![](joineRML_files/figure-html/re_comp_plot-1.png)
