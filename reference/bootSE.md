# Standard errors via bootstrap for an `mjoint` object

This function takes a model fit from an `mjoint` object and calculates
standard errors and confidence intervals for the main longitudinal and
survival coefficient parameters, including the latent association
parameters, using bootstrapping (Efron and Tibshirani, 2000).

## Usage

``` r
bootSE(
  object,
  nboot = 100,
  ci = 0.95,
  use.mle = TRUE,
  verbose = FALSE,
  control = list(),
  progress = TRUE,
  ncores = 1L,
  safe.boot = FALSE,
  ...
)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- nboot:

  the number of bootstrap samples. Default is `nboot=100`.

- ci:

  the confidence interval to be estimated using the percentile-method.
  Default is `ci=0.95` for a 95% confidence interval.

- use.mle:

  logical: should the algorithm use the maximizer from the converged
  model in `object` as initial values for coefficients in each bootstrap
  iteration. Default is `use.mle=TRUE`.

- verbose:

  logical: if `TRUE`, the parameter estimates and other convergence
  statistics are value are printed at each iteration of the MCEM
  algorithm. Default is `FALSE`.

- control:

  a list of control values with components:

  `nMC`

  :   integer: the initial number of Monte Carlo samples to be used for
      integration in the burn-in phase of the MCEM. Default is
      `nMC=`100*K*.

  `nMCscale`

  :   integer: the scale factor for the increase in Monte Carlo size
      when Monte Carlo has not reduced from the previous iteration.
      Default is `nMCscale=5`.

  `nMCmax`

  :   integer: the maximum number of Monte Carlo samples that the
      algorithm is allowed to reach. Default is `nMCmax=20000`.

  `burnin`

  :   integer: the number of iterations for 'burn-in' phase of the
      optimization algorithm. It is computationally inefficient to use a
      large number of Monte Carlo samples early on until one is
      approximately near the maximum likelihood estimate. Default is
      `burnin=`100*K* for `type='antithetic'` or `type='montecarlo'` and
      `burnin=`5 for `type='sobol'` or `type='halton'`. For standard
      methods, such a large burn-in will generally be unnecessary and
      can be reduced on an application-specific basis.

  `mcmaxIter`

  :   integer: the maximum number of MCEM algorithm iterations allowed.
      Default is `mcmaxIter=burnin+200`.

  `convCrit`

  :   character string: the convergence criterion to be used. See
      **Details**.

  `gammaOpt`

  :   character string: by default (`gammaOpt='NR'`), \\\gamma\\ is
      updated using a one-step Newton-Raphson iteration, with the
      Hessian matrix calculated exactly. If `gammaOpt='GN'`, a
      Gauss-Newton algorithm-type iteration is implemented, where the
      Hessian matrix is approximated based on calculations similar to
      those used for calculating the empirical information matrix? If it
      is used, then the step-length is adjusted by a nominal scaling
      parameter of 0.5 in order to reduce the chance of over-shooting
      the maximizer.

  `tol0`

  :   numeric: tolerance value for convergence in the parameters; see
      **Details**. Default is `tol0=1e-03`.

  `tol1`

  :   numeric: tolerance value for convergence in the parameters; see
      **Details**. Default is `tol1=1e-03`.

  `tol2`

  :   numeric: tolerance value for convergence in the parameters; see
      **Details**. Default is `tol2=5e-03` for `type='antithetic'` or
      `type='montecarlo'` and `tol2=1e-03` for `type='sobol'` or
      `type='halton'`.

  `tol.em`

  :   numeric: tolerance value for convergence in the multivariate
      linear mixed model (MV-LMM). When \\K \> 1\\, the optimal initial
      parameters are those from the MV-LMM, which is estimated using a
      separate EM algorithm. Since both the E- and M-steps are available
      in closed-form, this algorithm convergences relatively rapidly
      with a high precision. Default is min(`1e-04`, `tol2`).

  `rav`

  :   numeric: threshold when using `convCrit='sas'` that applies
      absolute change (when \\\<\\`rav`) or relative change (when
      \\\ge\\`rav`) criterion; see **Details**. Default is `0.1`, which
      is an order of magnitude higher than the SAS implementation.

  `type`

  :   character: type of Monte Carlo integration method to use. Options
      are

      `type='montecarlo'`

      :   Vanilla Monte Carlo sampling.

      `type='antithetic'`

      :   Variance reduction method using antithetic simulation. This is
          the default option.

      `type='sobol'`

      :   Quasi-Monte Carlo with a low deterministic Sobol sequence with
          Owen-type scrambling.

      `type='halton'`

      :   Quasi-Monte Carlo with a low deterministic Halton sequence.

- progress:

  logical: should a progress bar be shown on the console to indicate the
  percentage of bootstrap iterations completed? Default is
  `progress=TRUE`.

- ncores:

  integer: if more than one core is available, then the `bootSE`
  function can run in parallel via the
  [`foreach`](https://rdrr.io/pkg/foreach/man/foreach.html) function. By
  default, `ncores=1`, which defaults to serial mode. Note that if
  `ncores`\>1, then `progress` is set to `FALSE` by default, as it is
  not possible to display progress bars for parallel processes at the
  current time.

- safe.boot:

  logical: should each bootstrap replication be wrapped in a
  [`tryCatch`](https://rdrr.io/r/base/conditions.html) statement to
  catch errors (e.g. during the optimisation progress)? When model
  fitting throws errors, a new bootstrap sample is drawn for the current
  iteration and the model is re-fit; this process continues until a
  model fits successfully. Default is `FALSE`.

- ...:

  options passed to the `control` argument.

## Value

An object of class `bootSE`.

## Details

Standard errors and confidence intervals are obtained by repeated
fitting of the requisite joint model to bootstrap samples of the
original longitudinal and time-to-event data. Note that bootstrap is
done by sampling subjects, not individual records.

## Note

This function is computationally intensive. A dynamic progress bar is
displayed showing the percentage of bootstrap models fitted. On computer
systems with more than one core available, computational time can be
reduced by passing the argument `ncores` (with integer value \>1) to
`bootSE`, which implements parallel processing via the
[`foreach`](https://rdrr.io/pkg/foreach/man/foreach.html) package.
**Note:** if parallel processing is implemented, then the progress bar
is not displayed.

Due to random sampling, an `mjoint` model fitted to some bootstrap
samples may not converge within the specified control parameter
settings. The `bootSE` code discards any models that failed to converge
when calculating the standard errors and confidence intervals. If a
large proportion of models have failed to converge, it is likely that it
will need to be refitted with changes to the `control` arguments.

## References

Efron B, Tibshirani R. *An Introduction to the Bootstrap.* 2000; Boca
Raton, FL: Chapman & Hall/CRC.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
for approximate standard errors.

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a joint model with bivariate longitudinal outcomes

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    verbose = TRUE)

fit.boot <- bootSE(fit, 50, use.mle = TRUE, control = list(
    burnin = 25, convCrit = "either",
    tol0 = 6e-03, tol2 = 6e-03, mcmaxIter = 60))
} # }
```
