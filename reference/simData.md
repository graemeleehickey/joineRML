# Simulate data from a joint model

This function simulates multivariate longitudinal and time-to-event data
from a joint model.

## Usage

``` r
simData(
  n = 100,
  ntms = 5,
  beta = rbind(c(1, 1, 1, 1), c(1, 1, 1, 1)),
  gamma.x = c(1, 1),
  gamma.y = c(0.5, -1),
  sigma2 = c(1, 1),
  D = NULL,
  df = Inf,
  model = "intslope",
  theta0 = -3,
  theta1 = 1,
  censoring = TRUE,
  censlam = exp(-3),
  truncation = TRUE,
  trunctime = (ntms - 1) + 0.1
)
```

## Arguments

- n:

  the number of subjects to simulate data for.

- ntms:

  the maximum number of (discrete) time points to simulate repeated
  longitudinal measurements at.

- beta:

  a matrix of `dim=c(K,4)` specifying the coefficients of the fixed
  effects. The order in each row is intercept, time, a continuous
  covariate, and a binary covariate.

- gamma.x:

  a vector of `length=2` specifying the coefficients for the
  time-to-event baseline covariates, in the order of a continuous
  covariate and a binary covariate.

- gamma.y:

  a vector of `length=K` specifying the latent association parameters
  for each longitudinal outcome.

- sigma2:

  a vector of `length=K` specifying the residual standard errors.

- D:

  a positive-definite matrix specifying the variance-covariance matrix.
  If `model='int'`, the matrix has dimension `dim=c(K, K)`, else if
  `model='intslope'`, the matrix has dimension `dim =c(2K, 2K)`. If
  `D=NULL` (default), an identity matrix is assumed.

- df:

  a non-negative scalar specifying the degrees of freedom for the random
  effects if sampled from a multivariate *t*-distribution. The default
  is `df=Inf`, which corresponds to a multivariate normal distribution.

- model:

  follows the model definition in the
  [`joint`](https://rdrr.io/pkg/joineR/man/joint.html) function. See
  **Details** for choices.

- theta0, theta1:

  parameters controlling the failure rate. See Details.

- censoring:

  logical: if `TRUE`, includes an independent censoring time.

- censlam:

  a scale (\\\> 0\\) parameter for an exponential distribution used to
  simulate random censoring times for when `censoring=TRUE`.

- truncation:

  logical: if `TRUE`, adds a truncation time for a maximum event time.

- trunctime:

  a truncation time for use when `truncation=TRUE`.

## Value

A list of 2 `data.frame`s: one recording the requisite longitudinal
outcomes data, and one recording the time-to-event data.

## Details

The function `simData` simulates data from a joint model, similar to
that performed in Henderson et al. (2000). It works by first simulating
multivariate longitudinal data for all possible follow-up times using
random draws for the multivariate Gaussian random effects and residual
error terms. Data can be simulated assuming either random-intercepts
only in each of the longitudinal sub-models, or random-intercepts and
random-slopes. Currently, all models must have the same structure. The
failure times are simulated from proportional hazards time-to-event
models using the following methodologies:

- `model="int"`:

  The baseline hazard function is specified to be an exponential
  distribution with

  \$\$\lambda_0(t) = \exp{\theta_0}.\$\$

  Simulation is conditional on known time-independent effects, and the
  methodology of Bender et al. (2005) is used to simulate the failure
  time.

- `model="intslope"`:

  The baseline hazard function is specified to be a Gompertz
  distribution with

  \$\$\lambda_0(t) = \exp{\theta_0 + \theta_1 t}.\$\$

  In the usual representation of the Gompertz distribution, \\\theta_1\\
  is the shape parameter, and the scale parameter is equivalent to
  \\\exp(\theta_0)\\. Simulation is conditional on on a predictable
  (linear) time-varying process, and the methodology of Austin (2012) is
  used to simulate the failure time.

## References

Austin PC. Generating survival times to simulate Cox proportional
hazards models with time-varying covariates. *Stat Med.* 2012;
**31(29)**: 3946-3958.

Bender R, Augustin T, Blettner M. Generating survival times to simulate
Cox proportional hazards models. *Stat Med.* 2005; **24**: 1713-1723.

Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
measurements and event time data. *Biostatistics.* 2000; **1(4)**:
465-480.

## Author

Pete Philipson (<peter.philipson1@newcastle.ac.uk>) and Graeme L. Hickey
(<graemeleehickey@gmail.com>)

## Examples

``` r
beta <- rbind(c(0.5, 2, 1, 1),
c(2, 2, -0.5, -1))
D <- diag(4)
D[1, 1] <- D[3, 3] <- 0.5
D[1, 2] <- D[2, 1] <- D[3, 4] <- D[4, 3] <- 0.1
D[1, 3] <- D[3, 1] <- 0.01

sim <- simData(n = 250, beta = beta, D = D, sigma2 = c(0.25, 0.25),
               censlam = exp(-0.2), gamma.y = c(-.2, 1), ntms = 8)
#> 24.4% experienced event
```
