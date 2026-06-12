# Summary of an `mjoint` object

This function provides a summary of an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
summary(object, bootSE = NULL, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- bootSE:

  an object inheriting from class `bootSE` for the corresponding model.
  If `bootSE=NULL`, the function will attempt to utilize approximate
  standard error estimates (if available) calculated from the empirical
  information matrix.

- ...:

  additional arguments; currently none are used.

## Value

A list containing the coefficient matrices for the longitudinal and
time-to-event sub-models; variance-covariance matrix for the random
effects; residual error variances; log-likelihood of joint model; AIC
and BIC statistics; and model fit objects.

## References

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
measurements and event time data. *Biostatistics.* 2000; **1(4)**:
465-480.

Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the
joint analysis of time-to-event and multiple longitudinal variables.
*Stat Med.* 2002; **21**: 2369-2382.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`mjoint.object`](https://graemeleehickey.github.io/joineRML/reference/mjoint.object.md),
and [`summary`](https://rdrr.io/r/base/summary.html) for the generic
method description.

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
summary(fit2)
} # }
```
