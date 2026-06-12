# Extract an approximate variance-covariance matrix of estimated parameters from an `mjoint` object

Returns the variance-covariance matrix of the main parameters of a
fitted `mjoint` model object.

## Usage

``` r
# S3 method for class 'mjoint'
vcov(object, correlation = FALSE, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- correlation:

  logical: if `TRUE` returns the correlation matrix, otherwise returns
  the variance-covariance matrix (default).

- ...:

  additional arguments; currently none are used.

## Value

A variance-covariance matrix.

## Details

This is a generic function that extracts the variance-covariance matrix
of parameters from an `mjoint` model fit. It is based on a profile
likelihood, so no estimates are given for the baseline hazard function,
which is generally considered a nuisance parameter. It is based on the
empirical information matrix (see Lin et al. 2002, and McLachlan and
Krishnan 2008 for details), so is only approximate.

## Note

This function is not to be confused with
[`getVarCov`](https://rdrr.io/pkg/nlme/man/getVarCov.html), which
returns the extracted variance-covariance matrix for the random effects
distribution.

## References

Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the
joint analysis of time-to-event and multiple longitudinal variables.
*Stat Med.* 2002; **21**: 2369-2382.

McLachlan GJ, Krishnan T. *The EM Algorithm and Extensions*. Second
Edition. Wiley-Interscience; 2008.

## See also

[`vcov`](https://rdrr.io/r/stats/vcov.html) for the generic method
description, and [`cov2cor`](https://rdrr.io/r/stats/cor.html) for
details of efficient scaling of a covariance matrix into the
corresponding correlation matrix.

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

vcov(fit1)
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

vcov(fit2)
} # }
```
