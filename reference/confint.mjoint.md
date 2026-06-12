# Confidence intervals for model parameters of an `mjoint` object

This function computes confidence intervals for one or more parameters
in a fitted `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
confint(
  object,
  parm = c("Both", "Longitudinal", "Event"),
  level = 0.95,
  bootSE = NULL,
  ...
)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- parm:

  a character string specifying which sub-model parameter confidence
  intervals should be returned for. Can be specified as
  `parm='Longitudinal'` (multivariate longitudinal sub-model),
  `parm='Event'` (time-to-event sub-model), or `parm='both'` (default).

- level:

  the confidence level required. Default is `level=0.95` for a 95%
  confidence interval.

- bootSE:

  an object inheriting from class `bootSE` for the corresponding model.
  If `bootSE=NULL`, the function will attempt to utilize approximate
  standard error estimates (if available) calculated from the empirical
  information matrix.

- ...:

  additional arguments; currently none are used.

## Value

A matrix containing the confidence intervals for either the
longitudinal, time-to-event, or both sub-models.

## References

McLachlan GJ, Krishnan T. *The EM Algorithm and Extensions.* Second
Edition. Wiley-Interscience; 2008.

Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
measurements and event time data. *Biostatistics.* 2000; **1(4)**:
465-480.

Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the
joint analysis of time-to-event and multiple longitudinal variables.
*Stat Med.* 2002; **21**: 2369-2382.

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`bootSE`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md),
and [`confint`](https://rdrr.io/r/stats/confint.html) for the generic
method description.

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a classical univariate joint model with a single longitudinal outcome
# and a single time-to-event outcome

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

gamma <- c(0.1059417, 1.0843359)
sigma2 <- 0.03725999
beta <- c(4.9988669999, -0.0093527634, 0.0004317697)
D <- matrix(c(0.128219108, -0.006665505, -0.006665505, 0.002468688),
            nrow = 2, byrow = TRUE)

set.seed(1)
fit1 <- mjoint(formLongFixed = log.lvmi ~ time + age,
    formLongRandom = ~ time | num,
    formSurv = Surv(fuyrs, status) ~ age,
    data = hvd,
    timeVar = "time",
    inits = list(gamma = gamma, sigma2 = sigma2, beta = beta, D = D),
    control = list(nMCscale = 2, burnin = 5)) # controls for illustration only

confint(fit1, parm = "Longitudinal")
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
confint(fit2)
} # }
```
