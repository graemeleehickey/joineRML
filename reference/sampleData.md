# Sample from an `mjoint` object

Generic function used to sample a subset of data from an object of class
`mjoint` with a specific number of subjects.

## Usage

``` r
sampleData(object, size = NULL, replace = TRUE)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- size:

  number of subjects to include in the sampled subset. If `size=NULL`
  (default), then size is set equal to the number of subjects used to
  fit the `mjoint` model.

- replace:

  use replacement when sampling subjects? Default is `TRUE`. If
  replacement is used, then the subjects are re-labelled from 1 to
  `size`.

## Value

A list of 2 data.frames: one recording the requisite longitudinal
outcomes data, and the other recording the time-to-event data.

## Details

This function is primarily intended for internal use in the
[`bootSE`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
function in order to permit bootstrapping. However, it can be used for
other purposes given a fitted `mjoint` object.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md).

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
sampleData(fit2, size = 10)
} # }
```
