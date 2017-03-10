# joineRML 0.1.1.9000 (devel version).

## New features

* `bootSE()` now uses control parameters from fitted model and allows for individual parameter overwrite.

* Added a check that any initial covariance matrix given is positive-definite.

* Added a check that dimensions of any `inits` given match the calculated dimensions from the model formulae.

## Bug patches

* Patched a major bug in `gammaUpdate()` where ties in failure times were not being properly handled. The code for `gammaUpdate_approx()` always worked fine. This bug manifested when `bootSE()` was called due to the resampling with replacement yielding datasets with many more ties than in the original dataset used to fit the model. To fix it, the information matrix required scaling at each failure time by the number of failures in the data.

* Patched a minor bug with prevented univariate random-intercept models for being fitted.

* Patched a minor bug in plotting convergence traces.

* Patched a minor bug with bootstrapping univariate joint models without passing the MLEs as the initial values to the `mjoint()` call.

## Housekeeping

* Package now `Depends` on `survival` and `nlme` rather than `Imports` to allow `require()` statements to be removed from code.

* Prevented `roxygen` from exporting all functions. 

* Fixed Imports following CRAN Checks of v0.1.1.

* Minor documentation edits and corrections.

* Minor code tidy-up.

* LICENSE upgraded to GPL-3 to be compatible with `joineR` v1.1.0.

* Removed internal function `fast_nearPD()` from package code as was unused.

* Removed internal function `EexpArma()` from package code as was unused.

* Added unit tests.

# joineRML 0.1.1

* Patched sigma.R roxygen documentation to handle sigma S3 method changing from lme4 to stats in R v3.3.0.

# joineRML 0.1.0

* First package release.



