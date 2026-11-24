# Compare fitted joint models

Compare two or more fitted joint models using their joint
log-likelihoods, AIC, and BIC. For nested covariate models, also
calculate likelihood-ratio statistics and approximate chi-squared
p-values.

## Usage

``` r
# S3 method for class 'mjoint'
anova(object, ..., test = TRUE, nMC = 10000L, tol = 0.1)
```

## Arguments

- object:

  an object inheriting from class `mjoint`.

- ...:

  additional fitted `mjoint` objects, ordered from the reduced model to
  the full model when `test = TRUE`.

- test:

  logical: calculate likelihood-ratio tests between successive models?
  Default is `TRUE`. With `FALSE`, report the stored log-likelihoods,
  AIC, and BIC without testing or further integration.

- nMC:

  an integer of at least 2 giving the number of Halton integration
  points per subject used to re-evaluate each likelihood when testing.
  Default is `10000`. The first `floor(nMC / 2)` points are also used to
  assess numerical stability.

- tol:

  a positive finite tolerance for the absolute change in the
  likelihood-ratio statistic between the two integration sizes. Default
  is `0.1`. This is a numerical diagnostic, not a significance threshold
  or a bound on integration error.

## Value

An object of class `anova`, also inheriting from `data.frame`, with one
row per model. Columns `df`, `logLik`, `AIC`, and `BIC` describe each
model. When testing, `L.Ratio`, `Df`, and `Pr(>Chisq)` describe
comparisons with the preceding model, and `LR.change` is the absolute
change in the statistic between the two integration sizes. The
`integration` attribute records the method, integration sizes,
tolerance, and likelihoods at the smaller size.

## Details

All models must have converged and use the same subjects, longitudinal
responses and measurement times, and survival responses. Longitudinal
response expressions must also agree, so comparisons between different
response transformations are rejected. Models fitted with `pfs = FALSE`
are supported. Older fitted objects that lack the stored fitting
matrices must be refitted with the current version of `mjoint`.

Tests require unchanged random-effects and association structures, and
nested longitudinal fixed-effects and survival-covariate design
matrices. Each successive model must add at least one identifiable
coefficient. Nesting is checked using matrix column spaces, allowing
different factor contrasts. The reference distribution is chi-squared,
with degrees of freedom equal to the difference in parameter counts. The
common baseline hazard has the same event-time support in both models
and does not add to this difference. Tests that change the
random-effects structure, including tests of zero variance components,
are not supported. Non-nested models on the same observations can be
compared with `test = FALSE`.

When testing, the observed joint likelihood is re-evaluated at the
fitted coefficients using common deterministic Halton points for the
conditional Gaussian random effects. No models are refitted. The table's
log-likelihood, AIC, and BIC use these re-evaluated values, so they may
differ from the corresponding values stored in the fitted objects. BIC
uses the number of subjects. Neither the fitted objects nor R's
random-number state is changed.

A warning is issued and the p-value is set to `NA` if a likelihood-ratio
statistic is negative or changes by more than `tol` when the number of
integration points increases. Increase `nMC` to investigate integration
error. A negative statistic can also indicate that a fitted model needs
tighter convergence settings or better initial values. Passing this
diagnostic does not establish that either model reached the maximum, or
guarantee accurate integration or validity of the asymptotic
approximation.

A single-model call is not supported; no sequential term deletion or
automatic model fitting is performed.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`logLik.mjoint`](https://graemeleehickey.github.io/joineRML/reference/logLik.mjoint.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.lvmi), ]
reduced <- mjoint(
  formLongFixed = log.lvmi ~ time,
  formLongRandom = ~ 1 | num,
  formSurv = Surv(fuyrs, status) ~ age,
  data = hvd, timeVar = "time", seed = 123)
full <- update(reduced, formLongFixed = log.lvmi ~ time + sex)
anova(reduced, full)
# Check sensitivity to more integration points
anova(reduced, full, nMC = 20000)
# Descriptive comparison using the stored likelihoods
anova(reduced, full, test = FALSE)
} # }
```
