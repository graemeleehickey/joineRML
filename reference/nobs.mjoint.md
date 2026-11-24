# Extract the number of subjects

Extract the number of subjects

## Usage

``` r
# S3 method for class 'mjoint'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `mjoint` object.

- ...:

  Additional arguments (not used).

## Value

The number of subjects used to fit the joint model. Longitudinal
observation counts are reported separately by
[`glance()`](https://generics.r-lib.org/reference/glance.html).

## Details

This agrees with the `nobs` attribute of `logLik(object)` and the sample
size used in BIC.

## Author

Graeme L. Hickey
