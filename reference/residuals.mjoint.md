# Extract `mjoint` residuals

The residuals at level *i* are obtained by subtracting the fitted levels
at that level from the response vector.

## Usage

``` r
# S3 method for class 'mjoint'
residuals(object, level = 0, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- level:

  an optional integer giving the level of grouping to be used in
  extracting the residuals from object. Level values increase from
  outermost to innermost grouping, with level 0 corresponding to the
  population residuals and level 1 corresponding to subject-specific
  residuals. Defaults to `level=0`.

- ...:

  additional arguments; currently none are used.

## Value

A `list` of length *K* with each element a vector of residuals for the
*k*-th longitudinal outcome.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`fitted.mjoint`](https://graemeleehickey.github.io/joineRML/reference/fitted.mjoint.md)

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
