# Extract fixed effects estimates from an `mjoint` object

Extract fixed effects estimates from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
fixef(object, process = c("Longitudinal", "Event"), ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- process:

  character string: if `process='Longitudinal'` the fixed effects
  coefficients from the (multivariate) longitudinal sub-model are
  returned. Else, if `process='Event'`, the coefficients from the
  time-to-event sub-model are returned.

- ...:

  additional arguments; currently none are used.

## Value

A named vector of length equal to the number of sub-model coefficients
estimated.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

## See also

[`fixef`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) for the
generic method description, and
[`ranef.mjoint`](https://graemeleehickey.github.io/joineRML/reference/ranef.mjoint.md).

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

fixef(fit1, process = "Longitudinal")
fixef(fit1, process = "Event")
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

fixef(fit2, process = "Longitudinal")
fixef(fit2, process = "Event")
} # }
```
