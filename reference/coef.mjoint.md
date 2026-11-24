# Extract joint-model coefficients

Extract joint-model coefficients

## Usage

``` r
# S3 method for class 'mjoint'
coef(object, complete = TRUE, ..., process = NULL, k = NULL)
```

## Arguments

- object:

  A fitted `mjoint` object.

- complete:

  Logical; retained for compatibility with `coef.default()`.

- ...:

  Additional arguments (not supported).

- process:

  `NULL` (the default) returns the original coefficient list. Select
  `"Longitudinal"` or `"Event"` to return a named vector of regression
  coefficients.

- k:

  With `process = "Longitudinal"`, select one outcome. The default
  `NULL` returns coefficients for all longitudinal outcomes, retaining
  their numbered suffixes. Selecting an outcome removes those suffixes
  so names match its design matrix.

## Value

A list of all estimated parameters, or a named numeric vector for a
selected component.

## Details

The default return value is unchanged. Event coefficients include both
baseline covariates and latent association coefficients. The same
`process` and `k` selections can be supplied to
[`vcov()`](https://rdrr.io/r/stats/vcov.html) to obtain the
corresponding covariance matrix in the same order.

## Author

Graeme L. Hickey
