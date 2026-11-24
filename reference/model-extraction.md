# Extract model terms, frames, and design matrices

Extract model terms, frames, and design matrices

## Usage

``` r
# S3 method for class 'mjoint'
model.frame(formula, process = c("Longitudinal", "Event"), k = 1, ...)

# S3 method for class 'mjoint'
terms(x, process = c("Longitudinal", "Event"), k = 1, ...)

# S3 method for class 'mjoint'
model.matrix(object, process = c("Longitudinal", "Event"), k = 1, ...)
```

## Arguments

- process:

  Select `"Longitudinal"` (the default) or `"Event"`.

- k:

  The longitudinal outcome number, between 1 and the number of outcomes.
  Only used for `process = "Longitudinal"`.

- ...:

  Additional arguments. These methods extract the fitted model only;
  supplying other arguments, including new data, is an error.

- x, object, formula:

  A fitted `mjoint` object.

## Value

[`terms()`](https://rdrr.io/r/stats/terms.html) returns a `terms`
object; [`model.frame()`](https://rdrr.io/r/stats/model.frame.html)
returns a data frame;
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) returns a
numeric matrix.

## Details

The longitudinal methods describe the fixed-effects model for outcome
`k`. The event methods describe the baseline covariates in the survival
formula. Its design matrix has no intercept or latent association
columns. Rows follow the stored longitudinal subject-time order or
survival subject order.

Model frames contain the evaluated formula variables, including
transformed responses and covariates. Fitted terms and contrasts are
retained, so changing the original data or R's contrast options does not
change extraction. For objects fitted before version 0.4.9, frames are
reconstructed from stored data and terms when available.

Column names of a longitudinal design matrix match
`coef(object, process = "Longitudinal", k = k)`. Event design columns
match the baseline-covariate entries of
`coef(object, process = "Event")`; the latter also includes association
coefficients.

## Author

Graeme L. Hickey
