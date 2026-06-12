# Extract variance-covariance matrix of random effects from an `mjoint` object

Extract variance-covariance matrix of random effects from an `mjoint`
object.

## Usage

``` r
# S3 method for class 'mjoint'
getVarCov(obj, ...)
```

## Arguments

- obj:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- ...:

  additional arguments; currently none are used.

## Value

A variance-covariance matrix.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`getVarCov`](https://rdrr.io/pkg/nlme/man/getVarCov.html) for the
generic method description.

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

getVarCov(fit2)
} # }
```
