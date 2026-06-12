# Plot convergence time series for parameter vectors from an `mjoint` object

Plot convergence time series for parameter vectors from an `mjoint`
object.

## Usage

``` r
plotConvergence(object, params = "gamma", discard = FALSE)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- params:

  a string indicating what parameters are to be shown. Options are
  `params='gamma'` for the time-to-event sub-model covariate
  coefficients, including the latent association parameters;
  `params='beta'` for the longitudinal sub-model fixed effects
  coefficients; `params='sigma2'` for the residual error variances from
  the longitudinal sub-model; `params='D'` for the lower triangular
  matrix of the variance-covariance matrix of random effects;
  `params='loglik'` for the log-likelihood.

- discard:

  logical; if `TRUE` then the 'burn-in' phase iterations of the MCEM
  algorithm are discarded. Default is `discard=FALSE`.

## References

Wei GC, Tanner MA. A Monte Carlo implementation of the EM algorithm and
the poor man's data augmentation algorithms. *J Am Stat Assoc.* 1990;
**85(411)**: 699-704.

## See also

[`plot.mjoint`](https://graemeleehickey.github.io/joineRML/reference/plot.mjoint.md),
[`plot.default`](https://rdrr.io/r/graphics/plot.default.html),
[`par`](https://rdrr.io/r/graphics/par.html),
[`abline`](https://rdrr.io/r/graphics/abline.html).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)
