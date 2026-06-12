# Plot a `dynLong` object

Plots the conditional longitudinal expectations for a *new* subject
calculated using the
[`dynLong`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md)
function.

## Usage

``` r
# S3 method for class 'dynLong'
plot(x, main = NULL, xlab = NULL, ylab = NULL, grid = TRUE, estimator, ...)
```

## Arguments

- x:

  an object of class `dynLong` calculated by the
  [`dynLong`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md)
  function.

- main:

  an overall title for the plot: see
  [`title`](https://rdrr.io/r/graphics/title.html).

- xlab:

  a title for the x \[time\] axis: see
  [`title`](https://rdrr.io/r/graphics/title.html).

- ylab:

  a character vector of the titles for the *K* longitudinal outcomes
  y-axes: see [`title`](https://rdrr.io/r/graphics/title.html).

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

- ...:

  additional plotting arguments; currently limited to `lwd` and `cex`.
  See [`par`](https://rdrr.io/r/graphics/par.html) for details.

## Value

A dynamic prediction plot.

## References

Rizopoulos D. Dynamic predictions and prospective accuracy in joint
models for longitudinal and time-to-event data. *Biometrics*. 2011;
**67**: 819–829.

## See also

[`dynLong`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md)

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
out <- dynLong(fit2, hvd2)
plot(out, main = "Patient 1")
} # }
```
