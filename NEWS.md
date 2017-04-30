# joineRML 0.2.2

* This update is an attempt to overcome a FAIL status on the CRAN checks for macOS.

## Housekeeping

* Changed the technical vignette to Rnw (with engine = Sweave) from ltx (with engine = R.rsp) in attempt to remove some CMD check warnings.

* Shortened the vignette to make it compile quicker (removed execution of bootstrapping).

* Lots of tweaks of minor formatting tweaks in the documentation.

# joineRML 0.2.1

## Bug patches

* Fixed a bug that prevented factors with >2 levels being included in the time-to-event submodel.

* Fixed package registration, which strangely broke on R 3.4.0 for OSX platform.

## Housekeeping

* `joineRML` version 0.2.1 will depend on R version >=3.3.0 to remedy issue with `sigma.mjoint()` S3 issue.

* Added a new badge to the README.

# joineRML 0.2.0

## New features

* Add `residual()` and `fitted()` functions for `mjoint` objects.

* Added a plot function -- `plot.ranef.mjoint()` -- for `ranef.mjoint` objects.

* `bootSE()` now uses control parameters from fitted model and allows for individual parameter overwrite.

* Added a check that any initial covariance matrix given is positive-definite.

* Added a check that dimensions of any `inits` given match the calculated dimensions from the model formulae.

* Added a check that if multiple repeated longitudinal outcoems were given, that each subject contributes at least one measurement per outcome.

* Changed the API so that `postRE` and `approx.se` arguments are now replaced by a single `pfs` argument, which stands for post fit statistics. If `TRUE`, then *both* the approximate SEs and the BLUPs (and SEs) are calculated and returned. This change is to facilitate other post fit statistics, e.g. residuals.

* `sampleData()` now allows for sampling *without* replacement.

* `plot()` and `plotConvergence()` now have the option to discard burn-in phase iterations from the MCEM algorithm.

* `plot()` and `plotConvergence()` now plot the log-likelihood trace.

* `vcov()` now calculates the variance-covariance matrix rather than inside `mjoint()` and then extracting it. It also utilises the QR-decomposition inverse function and the Moore-Penrose matrix inverse, as in some cases the matrix was nearly singular.

* `hessian()` (and therefore `vcov()`) now calculate the contribution for the random effect variance terms rather than the random effect precision (1 divided by the variance) terms. The correct contribution to the score for off-diagonal terms is now also impleted.

* The left-hand side of the `formLongFixed` now handles transformations.

## Bug patches

* `vcov()` now returns the variance-covariance matrix as intended. Previously it was only returning the profile empirical information matrix.

* Patched a major bug in `gammaUpdate()` where ties in failure times were not being properly handled. The code for `gammaUpdate_approx()` always worked fine, as it was based only on the score vector. This bug manifested when `bootSE()` was called due to the resampling with replacement yielding datasets with many more ties than in the original dataset used to fit the model. To fix it, the information matrix required scaling at each failure time by the number of failures in the data. The formula for the information matrix in the Technical Details vignette has also been updated.

* Patched a minor bug with prevented univariate random-intercept models for being fitted.

* Patched a minor bug in plotting convergence traces.

* Patched a minor bug with bootstrapping univariate joint models without passing the MLEs as the initial values to the `mjoint()` call.

## Housekeeping

* Renamed `approxSE()` function to `hessian()`.

* Renamed `control` argument `earlyPhase` to `burnin`.

* Default settings for `control` parameters updated based on practical experience.

* Package now `Depends` on `survival` and `nlme` rather than `Imports` to allow `require()` statements to be removed from code.

* Prevented `roxygen` from exporting all functions. 

* Fixed Imports following CRAN Checks of v0.1.1.

* Minor documentation edits and corrections.

* Minor code tidy-up with slight speed-ups and stabilisations where found.

* LICENSE upgraded to GPL-3 to be compatible with `joineR` v1.1.0.

* Removed internal function `fast_nearPD()` from package code as was unused.

* Removed internal function `EexpArma()` from package code as was unused.

* Added unit tests.

* Registered native C++ routines and disabled symbol search to satisfy CRAN CMD checks.

# joineRML 0.1.1

* Patched sigma.R roxygen documentation to handle sigma S3 method changing from lme4 to stats in R v3.3.0.

# joineRML 0.1.0

* First package release.



