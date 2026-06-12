# Fitted `mjoint` object

An object returned by the `mjoint` function, inheriting from class
`mjoint` and representing a fitted joint model for multivariate
longitudinal and time-to-event data. Objects of this class have methods
for the generic functions `coef`, `logLik`, `plot`, `print`, `ranef`,
`fixef`, `summary`, `AIC`, `getVarCov`, `vcov`, `confint`, `sigma`,
`fitted`, `residuals`, and `formula`.

## Usage

``` r
mjoint.object
```

## Format

An object of class `NULL` of length 0.

## Value

A list with the following components.

- `coefficients`:

  a list with the estimated coefficients. The components of this list
  are:

  `beta`

  :   the vector of fixed effects for the linear mixed effects
      sub-model.

  `D`

  :   the variance-covariance matrix of the random effects.

  `sigma2`

  :   the measurement error standard deviations for the linear mixed
      effects sub-model.

  `haz`

  :   the estimated baseline hazard values for each unique failure time.
      Note that this is the *centered* hazard, equivalent to that
      returned by
      [`coxph.detail`](https://rdrr.io/pkg/survival/man/coxph.detail.html).

  `gamma`

  :   the vector of baseline covariates for the survival model and the
      latent association coefficient parameter estimates.

- `history`:

  a matrix with parameter estimates at each iteration of the MCEM
  algorithm.

- `nMC.hx`:

  a vector with the number of Monte Carlo samples for each MCEM
  algorithm iteration.

- `formLongFixed`:

  a list of formulae for the fixed effects component of each
  longitudinal outcome.

- `formLongRandom`:

  a list of formulae for the fixed effects component of each
  longitudinal outcome. The length of the list will be equal to
  `formLongFixed`.

- `formSurv`:

  a formula specifying the proportional hazards regression model (not
  including the latent association structure).

- `data`:

  a list of data.frames for each longitudinal outcome.

- `survData`:

  a data.frame of the time-to-event dataset.

- `timeVar`:

  a character string vector of length K denoting the column name(s) for
  time in `data`.

- `id`:

  a character string denoting the column name for subject IDs in `data`
  and `survData`.

- `dims`:

  a list giving the dimensions of model parameters with components:

  `p`

  :   a vector of the number of fixed effects for each longitudinal
      outcome.

  `r`

  :   a vector of the number of random effects for each longitudinal
      outcome.

  `K`

  :   an integer of the number of different longitudinal outcome types.

  `q`

  :   an integer of the number of baseline covariates in the
      time-to-event sub-model.

  `n`

  :   an integer of the total number of subjects in the study.

  `nk`

  :   a vector of the number of measurements for each longitudinal
      outcome.

- `sfit`:

  an object of class `coxph` for the separate time-to-event model fit.
  See [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for
  details.

- `lfit`:

  a list of objects each of class `lme` from fitting separate linear
  mixed effects models; one per each longitudinal outcome type. See
  [`lme`](https://rdrr.io/pkg/nlme/man/lme.html) for details.

- `log.lik0`:

  the combined log-likelihood from separate sub-model fits.

- `log.lik`:

  the log-likelihood from the joint model fit.

- `ll.hx`:

  a vector of the log-likelihood values for each MCEM algorithm
  interaction.

- `control`:

  a list of control parameters used in the estimation of the joint
  model. See
  [`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  for details.

- `finalnMC`:

  the final number of Monte Carlo samples required prior to convergence.

- `call`:

  the matched call.

- `inits`:

  the initial values passed as an argument in the `mjoint` function.

- `inits.long`:

  the computed initial values from fitting a multivariate longitudinal
  model.

- `inits.surv`:

  the computed initial values from fitting a Cox proportional hazards
  model with time-dependent covariates calculated from the fitted
  multivariate LME model.

- `conv`:

  logical: did the MCEM algorithm converge within the specified maximum
  number of iterations?

- `comp.time`:

  a vector of length 2 with each element an object of class `difftime`
  that reports the *total* time taken for model fitting (including all
  stages) and the time spent in the *EM algorithm*.

## Post model fit statistics

If `pfs=TRUE`, indicating that post-fit statistics are to be returned,
then the output also includes the following objects.

- `vcov`:

  the variance-covariance matrix of model parameters, as approximated by
  the empirical information matrix, is reported. See
  [`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  for details.

- `SE.approx`:

  the square-root of the diagonal of `vcov` is returned, which are
  estimates of the standard errors for the parameters.

- `Eb`:

  a matrix with the estimated random effects values for each subject.

- `Vb`:

  an array with the estimated variance-covariance matrices for the
  random effects values for each subject.

- `dmats`:

  a list of length 3 containing the design matrices, data frames, and
  vectors used in the MCEM algorithm. These are required for prediction
  and to calculate the residuals and . The 3 items in the list are `l`
  (longitudinal data), `t` (time-to-event data), and `z` (design
  matrices expanded over unique failure times). These are not intended
  to be extracted by the user.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
