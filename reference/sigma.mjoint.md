# Extract residual standard deviation(s) from an `mjoint` object

Extract residual standard deviation(s) from an `mjoint` object.

## Usage

``` r
# S3 method for class 'mjoint'
sigma(object, ...)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- ...:

  additional arguments; currently none are used.

## Value

a number (standard deviation) if \\K = 1\\ (univariate model), or a
vector if \\K\>1\\ (multivariate model).

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`sigma`](https://rdrr.io/pkg/lme4/man/sigma.html) in the **lme4**
package.

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
