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

## Details

These values are the square roots of the residual variances stored in
`coef(object)$sigma2`; they are not standard errors of estimated
parameters. Names identify the longitudinal outcomes, or use `sigma_1`,
`sigma_2`, and so on when the formula list is unnamed.

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`sigma`](https://rdrr.io/pkg/lme4/man/sigma.html) in the **lme4**
package.

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
