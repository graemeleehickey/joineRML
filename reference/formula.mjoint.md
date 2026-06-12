# Extract model formulae from an `mjoint` object

Extract model formulae from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
formula(x, process = c("Longitudinal", "Event"), k = 1, ...)
```

## Arguments

- x:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- process:

  character string: if `process='Longitudinal'` a fixed effects formula
  from the (multivariate) longitudinal sub-model is returned for the
  `k`-th outcome. Else, if `process='Event'`, the time-to-event model
  formula is returned.

- k:

  integer: a number between 1 and *K* (the total number of longitudinal
  outcomes) that specifies the longitudinal outcome of interest.

- ...:

  additional arguments; currently none are used.

## Value

An object of class "formula" which contains a symbolic model formula for
the separate sub-model fixed effect terms only.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.

## See also

[`formula`](https://rdrr.io/r/stats/formula.html) for the generic method
description, and
[`ranef.mjoint`](https://graemeleehickey.github.io/joineRML/reference/ranef.mjoint.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
