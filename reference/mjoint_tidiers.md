# Tidying methods for joint models for time-to-event data and multivariate longitudinal data

These methods tidy the coefficients of joint models for time-to-event
data and multivariate longitudinal data of the `mjoint` class from the
`joineRML` package.

## Usage

``` r
# S3 method for class 'mjoint'
tidy(
  x,
  component = "survival",
  bootSE = NULL,
  conf.int = FALSE,
  conf.level = 0.95,
  ...
)

# S3 method for class 'mjoint'
augment(x, data = x$data, ...)

# S3 method for class 'mjoint'
glance(x, ...)
```

## Arguments

- x:

  An object of class `mjoint`.

- component:

  Either `survival` (the survival component of the model, default) or
  `longitudinal` (the longitudinal component).

- bootSE:

  An object of class `bootSE` for the corresponding model. If
  `bootSE = NULL` (the default), the function will use approximate
  standard error estimates calculated from the empirical information
  matrix.

- conf.int:

  Include (1 - `conf.level`)% confidence intervals? Defaults to `FALSE`.

- conf.level:

  The confidence level required.

- ...:

  extra arguments (not used)

- data:

  Original data this was fitted on, in a list (e.g. `list(data)`). This
  will be extracted from `x` if not given.

## Value

All tidying methods return a `data.frame` without rownames. The
structure depends on the method chosen.

`tidy` returns one row for each estimated fixed effect depending on the
`component` parameter. It contains the following columns:

- term:

  The term being estimated

- estimate:

  Estimated value

- std.error:

  Standard error

- statistic:

  Z-statistic

- p.value:

  P-value computed from Z-statistic

- conf.low:

  The lower bound of a confidence interval on `estimate`, if required

- conf.high:

  The upper bound of a confidence interval on `estimate`, if required

.

`augment` returns one row for each original observation, with columns
(each prepended by a .) added. Included are the columns:

- .fitted_j_0:

  population-level fitted values for the j-th longitudinal process

- .fitted_j_1:

  individuals-level fitted values for the j-th longitudinal process

- .resid_j_0:

  population-level residuals for the j-th longitudinal process

- .resid_j_1:

  individual-level residuals for the j-th longitudinal process

See
[`fitted.mjoint`](https://graemeleehickey.github.io/joineRML/reference/fitted.mjoint.md)
and
[`residuals.mjoint`](https://graemeleehickey.github.io/joineRML/reference/residuals.mjoint.md)
for more information on the difference between population-level and
individual-level fitted values and residuals.

`glance` returns one row with the columns

- sigma2_j:

  the square root of the estimated residual variance for the j-th
  longitudinal process

- AIC:

  the Akaike Information Criterion

- BIC:

  the Bayesian Information Criterion

- logLik:

  the data's log-likelihood under the model

.

## Note

If fitting a joint model with a single longitudinal process, please make
sure you are using a named `list` to define the formula for the fixed
and random effects of the longitudinal submodel.

## Author

Alessandro Gasparini (<alessandro.gasparini@ki.se>)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a joint model with bivariate longitudinal outcomes
library(joineRML)
data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) &
  !is.na(heart.valve$log.lvmi) &
  heart.valve$num <= 50, ]
fit <- mjoint(
  formLongFixed = list(
    "grad" = log.grad ~ time + sex + hs,
    "lvmi" = log.lvmi ~ time + sex
  ),
  formLongRandom = list(
    "grad" = ~ 1 | num,
    "lvmi" = ~ time | num
  ),
  formSurv = Surv(fuyrs, status) ~ age,
  data = hvd,
  inits = list("gamma" = c(0.11, 1.51, 0.80)),
  timeVar = "time"
)

# Extract the survival fixed effects
tidy(fit)

# Extract the longitudinal fixed effects
tidy(fit, component = "longitudinal")

# Extract the survival fixed effects with confidence intervals
tidy(fit, ci = TRUE)

# Extract the survival fixed effects with confidence intervals based on 
# bootstrapped standard errors
bSE <- bootSE(fit, nboot = 5, safe.boot = TRUE)
tidy(fit, bootSE = bSE, ci = TRUE)

# Augment original data with fitted longitudinal values and residuals
hvd2 <- augment(fit)

# Extract model statistics
glance(fit)
} # }
```
