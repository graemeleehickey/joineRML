# Plot a `dynSurv` object

Plots the conditional time-to-event distribution for a *new* subject
calculated using the
[`dynSurv`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
function.

## Usage

``` r
# S3 method for class 'dynSurv'
plot(
  x,
  main = NULL,
  xlab = NULL,
  ylab1 = NULL,
  ylab2 = NULL,
  grid = TRUE,
  estimator,
  smooth = FALSE,
  ...
)
```

## Arguments

- x:

  an object of class `dynSurv` calculated by the
  [`dynSurv`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
  function.

- main:

  an overall title for the plot: see
  [`title`](https://rdrr.io/r/graphics/title.html).

- xlab:

  a title for the x \[time\] axis: see
  [`title`](https://rdrr.io/r/graphics/title.html).

- ylab1:

  a character vector of the titles for the *K* longitudinal outcomes
  y-axes: see [`title`](https://rdrr.io/r/graphics/title.html).

- ylab2:

  a title for the event-time outcome axis: see
  [`title`](https://rdrr.io/r/graphics/title.html).

- grid:

  adds a rectangular grid to an existing plot: see
  [`grid`](https://rdrr.io/r/graphics/grid.html).

- estimator:

  a character string that can take values `mean` or `median` to specify
  what prediction statistic is plotted from an objecting inheritting of
  class `dynSurv`. Default is `estimator='median'`. This argument is
  ignored for non-simulated `dynSurv` objects, i.e. those of
  `type='first-order'`, as in that case a mode-based prediction is
  plotted.

- smooth:

  logical: whether to overlay a smooth survival curve (see **Details**).
  Defaults to `FALSE`.

- ...:

  additional plotting arguments; currently limited to `lwd` and `cex`.
  See [`par`](https://rdrr.io/r/graphics/par.html) for details.

## Value

A dynamic prediction plot.

## Details

The `joineRML` package is based on a semi-parametric model, such that
the baseline hazards function is left unspecified. For prediction, it
might be preferable to have a smooth survival curve. Rather than
changing modelling framework *a prior*, a constrained B-splines
non-parametric median quantile curve is estimated using
[`cobs`](https://rdrr.io/pkg/cobs/man/cobs.html), with a penalty
function of \\\lambda=1\\, and subject to constraints of monotonicity
and \\S(t)=1\\.

## References

Ng P, Maechler M. A fast and efficient implementation of qualitatively
constrained quantile smoothing splines. *Statistical Modelling*. 2007;
**7(4)**: 315-328.

Rizopoulos D. Dynamic predictions and prospective accuracy in joint
models for longitudinal and time-to-event data. *Biometrics*. 2011;
**67**: 819–829.

## See also

[`dynSurv`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)

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

hvd2 <- droplevels(hvd[hvd$num == 1, ])
out1 <- dynSurv(fit2, hvd2)
plot(out1, main = "Patient 1")
} # }

if (FALSE) { # \dontrun{
# Monte Carlo simulation with 95% confidence intervals on plot

out2 <- dynSurv(fit2, hvd2, type = "simulated", M = 200)
plot(out2, main = "Patient 1")
} # }
```
