# Extract log-likelihood from an `mjoint` object

Extract log-likelihood from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
logLik(object, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- ...:

  additional arguments; currently none are used.

## Value

Returns an object of class `logLik`. This is a number with two
attributes: `df` (degrees of freedom), giving the number of parameters
in the model, and `nobs`, the number of observations used in estimation.

## References

Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
measurements and event time data. *Biostatistics.* 2000; **1(4)**:
465-480.

## See also

[`logLik`](https://rdrr.io/r/stats/logLik.html) for the generic method
description.

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

logLik(fit2)
} # }
```
