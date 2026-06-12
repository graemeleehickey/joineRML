# Plot diagnostics from an `mjoint` object

Plot diagnostics from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
plot(x, type = "convergence", ...)
```

## Arguments

- x:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- type:

  currently the only option is `type='convergence'` for graphical
  examination of convergence over MCEM iteration.

- ...:

  other parameters passed to
  [`plotConvergence`](https://graemeleehickey.github.io/joineRML/reference/plotConvergence.md).

## See also

[`plot.default`](https://rdrr.io/r/graphics/plot.default.html),
[`par`](https://rdrr.io/r/graphics/par.html),
[`abline`](https://rdrr.io/r/graphics/abline.html).

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

plot(fit1, param = "beta")  # LMM fixed effect parameters
plot(fit1, param = "gamma") # event model parameters
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
    control = list(burnin = 50),
    verbose = TRUE)

plot(fit2, type = "convergence", params = "gamma")
} # }
```
