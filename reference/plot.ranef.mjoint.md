# Plot a `ranef.mjoint` object

Displays a plot of the BLUPs and approximate 95% prediction interval for
each subject.

## Usage

``` r
# S3 method for class 'ranef.mjoint'
plot(x, ...)
```

## Arguments

- x:

  an object inheriting from class `ranef.mjoint`, representing the
  estimated random effects for the `mjoint` object from which it was
  produced.

- ...:

  additional arguments; currently none are used.

## Value

an object inheriting from class `ggplot`, which displays a trellis plot
with a separate panel for each effect, showing a dotplot (with optional
error bars indicating approximate 95% prediction intervals if the
argument `postVar=TRUE` is set in the call to
[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html)) for each
subject (by row).

## References

Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS.* New York:
Springer Verlag; 2000.

## See also

[`ranef.mjoint`](https://graemeleehickey.github.io/joineRML/reference/ranef.mjoint.md).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r
if (FALSE) { # \dontrun{
require(ggplot2)
data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]
set.seed(1)

fit1 <- mjoint(formLongFixed = log.lvmi ~ time,
    formLongRandom = ~ time | num,
    formSurv = Surv(fuyrs, status) ~ 1,
    data = hvd,
    timeVar = "time")

plot(ranef(fit1, postVar = TRUE))
} # }
```
