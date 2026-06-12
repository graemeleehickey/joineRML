# Fit a joint model to time-to-event data and multivariate longitudinal data

This function fits the joint model proposed by Henderson et al. (2000),
but extended to the case of multiple continuous longitudinal measures.
The time-to-event data is modelled using a Cox proportional hazards
regression model with time-varying covariates. The multiple longitudinal
outcomes are modelled using a multivariate version of the Laird and Ware
linear mixed model. The association is captured by a multivariate latent
Gaussian process. The model is estimated using a (quasi-) Monte Carlo
Expectation Maximization (MCEM) algorithm.

## Usage

``` r
mjoint(
  formLongFixed,
  formLongRandom,
  formSurv,
  data,
  survData = NULL,
  timeVar,
  inits = NULL,
  verbose = FALSE,
  pfs = TRUE,
  control = list(),
  ...
)
```

## Arguments

- formLongFixed:

  a list of formulae for the fixed effects component of each
  longitudinal outcome. The left-hand side defines the response, and the
  right-hand side specifies the fixed effect terms. If a single formula
  is given (either as a list of length 1 or a formula), then it is
  assumed that a standard univariate joint model is being fitted.

- formLongRandom:

  a list of one-sided formulae specifying the model for the random
  effects effects of each longitudinal outcome. The length of the list
  must be equal to `formLongFixed`.

- formSurv:

  a formula specifying the proportional hazards regression model (not
  including the latent association structure). See
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for examples.

- data:

  a list of `data.frame` objects for each longitudinal outcome in which
  to interpret the variables named in the `formLongFixed` and
  `formLongRandom`. The `list` structure enables one to include multiple
  longitudinal outcomes with different measurement protocols. If the
  multiple longitudinal outcomes are measured at the same time points
  for each patient, then a `data.frame` object can be given instead of a
  `list`. It is assumed that each data frame is in long format. `tibble`
  objects are automatically converted to plain `data.frame` objects.

- survData:

  a `data.frame` in which to interpret the variables named in the
  `formSurv`. This is optional, and if not given, the required data is
  searched for in `data`. Default is `survData=NULL`.

- timeVar:

  a character string indicating the time variable in the linear mixed
  effects model. If there are multiple longitudinal outcomes and the
  time variable is labelled differently in each model, then a character
  string vector can be given instead.

- inits:

  a list of initial values for some or all of the parameters estimated
  in the model. Default is `NULL`, with initial values estimated using
  separate multivariate linear mixed effects and Cox proportional hazard
  regression models.

- verbose:

  logical: if `TRUE`, the parameter estimates and other convergence
  statistics are value are printed at each iteration of the MCEM
  algorithm. Default is `FALSE`.

- pfs:

  logical: if `TRUE`, then assuming the MCEM algorithm has converged,
  post-fit statistics including the posterior means and variances of the
  random effects, and the approximate standard errors are calculated and
  returned as part of the model object. Default is `TRUE`. If `FALSE`,
  then these additional calculations are not performed, which can reduce
  the overall computational time. This option is intended to be used
  with computationally intensive routines such as simulation and
  bootstrap standard error estimation where these calculations are not
  required.

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

- ...:

  options passed to the `control` argument.

## Value

An object of class `mjoint`. See
[`mjoint.object`](https://graemeleehickey.github.io/joineRML/reference/mjoint.object.md)
for details.

## Details

Function `mjoint` fits joint models for time-to-event and multivariate
longitudinal data. A review of relevant statistical methodology for
joint models of multivariate data is given in Hickey et al. (2016). This
is a direct extension of the models developed in the seminal works of
Wulfsohn and Tsiatis (1997) and Henderson et al. (2000), with the
calculations based largely on Lin et al. (2002) who also extended the
model to multivariate joint data. A more detailed explanation of the
model formulation is given in the technical vignette. Each longitudinal
outcome is modelled according to a linear mixed model (LMM), akin to the
models fit by [`lme`](https://rdrr.io/pkg/nlme/man/lme.html), with
independent and identically distributed Gaussian errors. The latent term
in each model (specified by `formLongRandom`) is a linear combination of
subject-specific zero-mean Gaussian random effects with outcome-specific
variance components. We denote these as \\W\_{i1}(t, b\_{i1})\\,
\\W\_{i2}(t, b\_{i2})\\, \\\dots\\, \\W\_{iK}(t, b\_{iK})\\, for the
*K*-outcomes. Usually, \\W(t, b)\\ will be specified as either \\b_0\\
(a random-intercepts model) or \\b_0 + b_1t\\ (a random-intercepts and
random-slopes model); however, more general structures are allowed The
time-to-event model is modelled as per the usual Cox model formulation,
with an additional (possibly) time-varying term given by

\$\$\gamma\_{y1} W\_{i1}(t, b\_{i1}) + \gamma\_{y2} W\_{i2}(t,
b\_{i2}) + \dots + \gamma\_{yK} W\_{iK}(t, b\_{iK}),\$\$

where \\\gamma_y\\ is a parameter vector of proportional latent
association parameters of length *K* for estimation.

The optimization routine is based on a Monte Carlo Expectation
Maximization algorithm (MCEM) algorithm, as described by Wei and Tanner
(1990). With options for using antithetic simulation for variance
reduction in the Monte Carlo integration, or quasi-Monte Carlo based on
low order deterministic sequences.

## Convergence criteria

The routine internally scales and centers data to avoid overflow in the
argument to the exponential function. These actions do not change the
result, but lead to more numerical stability. Several convergence
criteria are available:

- `abs`:

  the maximum absolute parameter change is \\\<\\`tol0`. The baseline
  hazard parameters are not included in this convergence statistic.

- `rel`:

  the maximum (absolute) relative parameter change is \\\<\\`tol2`. A
  small value (`tol1`) is added to the denominator of the relative
  change statistic to avoid numerical problems when the parameters are
  close to zero.

- `either`:

  *either* the `abs` or `rel` criteria are satisfied.

- `sas`:

  if \\\| \theta_p \| \< \\`rav`, then the `abs` criteria is applied for
  the *l*-th parameter; otherwise, `rel` is applied. This is the
  approach used in the SAS EM algorithm program for missing data.

Due to the Monte Caro error, the algorithm could spuriously declare
convergence. Therefore, we require convergence to be satisfied for 3
consecutive iterations. The algorithm starts with a low number of Monte
Carlo samples in the 'burn-in' phase, as it would be computationally
inefficient to use a large sample whilst far away from the true
maximizer. After the algorithm moves out of this adaptive phase, it uses
an automated criterion based on the coefficient of variation of the
relative parameter change of the last 3 iterations to decide whether to
increase the Monte Carlo sample size. See the technical vignette and
Ripatti et al. (2002) for further details.

## Standard error estimation

Approximate standard errors (SEs) are calculated (if `pfs=TRUE`). These
are based on the empirical observed information function (McLachlan &
Krishnan, 2008). Through simulation studies, we have found that this
approximation does not work particularly well for \\n \< 100\\ (where
*n* is the number of subjects). In these cases, one would need to appeal
to the bootstrap SE estimation approach. However, in practice, the
reliability of the approximate SEs will depend of a multitude of
factors, including but not limited to, the average number of repeated
measurements per subject, the total number of events, and the
convergence of the MCEM algorithm.

Bootstrap SEs are also available, however they are not calculated using
the `mjoint` function due to the intense computational time. Instead, a
separate function is available: `bootSE`, which takes the fitted joint
model as its main argument. Given a fitted joint model (of class
`mjoint`) and a bootstrap fit object (of class `bootSE`), the SEs
reported in the model can be updated by running
`summary(fit_obj, bootSE = boot_obj)`. For details, consult the
[`bootSE`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
documentation.

## References

Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
measurements and event time data. *Biostatistics.* 2000; **1(4)**:
465-480.

Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. Joint modelling
of time-to-event and multivariate longitudinal outcomes: recent
developments and issues. *BMC Med Res Methodol.* 2016; **16(1)**: 117.

Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the
joint analysis of time-to-event and multiple longitudinal variables.
*Stat Med.* 2002; **21**: 2369-2382.

McLachlan GJ, Krishnan T. *The EM Algorithm and Extensions.* Second
Edition. Wiley-Interscience; 2008.

Ripatti S, Larsen K, Palmgren J. Maximum likelihood inference for
multivariate frailty models using an automated Monte Carlo EM algorithm.
*Lifetime Data Anal.* 2002; **8**: 349-360.

Wei GC, Tanner MA. A Monte Carlo implementation of the EM algorithm and
the poor man's data augmentation algorithms. *J Am Stat Assoc.* 1990;
**85(411)**: 699-704.

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

## See also

[`mjoint.object`](https://graemeleehickey.github.io/joineRML/reference/mjoint.object.md),
[`bootSE`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md),
[`plot.mjoint`](https://graemeleehickey.github.io/joineRML/reference/plot.mjoint.md),
[`summary.mjoint`](https://graemeleehickey.github.io/joineRML/reference/summary.mjoint.md),
[`getVarCov.mjoint`](https://graemeleehickey.github.io/joineRML/reference/getVarCov.mjoint.md),
[`simData`](https://graemeleehickey.github.io/joineRML/reference/simData.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a classical univariate joint model with a single longitudinal outcome
# and a single time-to-event outcome

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

set.seed(1)
fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
    formLongRandom = ~ time | num,
    formSurv = Surv(fuyrs, status) ~ age,
    data = hvd,
    timeVar = "time",
    control = list(nMCscale = 2, burnin = 5)) # controls for illustration only
summary(fit1)
} # }

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
summary(fit2)
} # }

if (FALSE) { # \dontrun{
# Fit a univariate joint model and compare to the joineR package

data(pbc2)
pbc2$log.b <- log(pbc2$serBilir)

# joineRML package
fit.joineRML <- mjoint(
    formLongFixed = list("log.bil" = log.b ~ year),
    formLongRandom = list("log.bil" = ~ 1 | id),
    formSurv = Surv(years, status2) ~ age,
    data = pbc2,
    timeVar = "year")
summary(fit.joineRML)

# joineR package
pbc.surv <- UniqueVariables(pbc2, var.col = c("years", "status2"),
                            id.col = "id")
pbc.long <- pbc2[, c("id", "year", "log.b")]
pbc.cov <- UniqueVariables(pbc2, "age", id.col = "id")
pbc.jd <- jointdata(longitudinal = pbc.long, baseline = pbc.cov,
                    survival = pbc.surv, id.col = "id", time.col = "year")
fit.joineR <- joint(data = pbc.jd,
    long.formula = log.b ~ 1 + year,
    surv.formula = Surv(years, status2) ~ age,
    model = "int")
summary(fit.joineR)
} # }
```
