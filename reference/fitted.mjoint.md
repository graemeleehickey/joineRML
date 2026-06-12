# Extract `mjoint` fitted values

The fitted values at level *i* are obtained by adding together the
population fitted values (based only on the fixed effects estimates) and
the estimated contributions of the random effects.

## Usage

``` r
# S3 method for class 'mjoint'
fitted(object, level = 0, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- level:

  an optional integer giving the level of grouping to be used in
  extracting the fitted values from object. Level values increase from
  outermost to innermost grouping, with level 0 corresponding to the
  population fitted values and level 1 corresponding to subject-specific
  fitted values Defaults to level 0.

- ...:

  additional arguments; currently none are used.

## Value

A `list` of length *K* with each element a vector of fitted values for
the *k*-th longitudinal outcome.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
[`residuals.mjoint`](https://graemeleehickey.github.io/joineRML/reference/residuals.mjoint.md)

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
