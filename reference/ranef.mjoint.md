# Extract random effects estimates from an `mjoint` object

Extract random effects estimates from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
ranef(object, postVar = FALSE, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- postVar:

  logical: if `TRUE` the variance of the posterior distribution is also
  returned.

- ...:

  additional arguments; currently none are used.

## Value

A `data.frame` (also of class `ranef.mjoint`) with rows denoting the
individuals and columns the random effects (e.g., intercepts, slopes,
etc.). If `postVar=TRUE`, the numeric matrix has an extra attribute,
`postVar`.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

## See also

[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html) for the
generic method description, and
[`fixef.mjoint`](https://graemeleehickey.github.io/joineRML/reference/fixef.mjoint.md).
To plot `ranef.mjoint` objects, see
[`plot.ranef.mjoint`](https://graemeleehickey.github.io/joineRML/reference/plot.ranef.mjoint.md).

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

ranef(fit2)
} # }
```
