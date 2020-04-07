# joineRML 0.4.4

## New features

* Added `broom` compatability.

## Bug patches

* Fixed issue that required subject IDs to be first column in dataset.

## Housekeeping

* Added @ellessenne as package author.
* Updates to DESCRIPTION.
* Documentation updates.

# joineRML 0.4.3

## Housekeeping

* Graeme Hickey taken over the package as creator and maintainer.
* Updates to DESCRIPTION
* Updates to Makevars and Makevars.win
* Updates to appveyor.yml and .travis.yml settings
* Updates to .gitignore

# joineRML 0.4.2

## Bug patches

* Fixed issue with citation date due to DATE not being present in the DESCRIPTION file.

## Housekeeping

* Added Zenodo DOI badge to README.

* Added ORCID IDs for authors.

* Remove unused objects once finished with them to try and reduce memory overheads.

* Changed maintainer to Pete Philipson.

# joineRML 0.4.1

## New features

* Added smoothed predicted survival curves to the `plot.dynSuv()`. Smoothing is based on the constrained B-splines method.

* `dynSurv()` now includes an argument to specify a horizon time from the last known observation time.

* `simData()` now includes an argument to choose multivariate *t*-distributed random effects with varying degrees of freedom, thus allowing for sensitivity analyses of heavier tail distributions.

## Housekeeping

* Minor corrections to documentation.

* Minor bug fixes to plotting functions.

* Added MRC to DESCRIPTION as funder.

* Added Depsy badge to README.

# joineRML 0.4.0

## New features

* The MCEM algorithm can be implemented with vanilla Monte Carlo and quasi-Monte Carlo (using either the scrambled Sobol sequence or the Halton sequence). This is implemented through the `type` argument passed to the list of `control` parameters in `mjoint()`.

* `bootSE()` now has option to use parallel computing via the `foreach` package.

## Bug patches

* Fixed some small errors in the `epileptic.qol` dataset.

* Fixed situation where a tibble might be given as the dataset (#55 @ellessenne).

* Catches errors in bootstrap due to "bad" data and automatically restarts the bootstrap (#57 @ellessenne) 

## Housekeeping

* Added hex sticker badge.

* Moved the make files and raw data for `qol` and `renal` datasets into the `~/data-raw/` directory.

* Added @ellessenne as package contributor for testing and bug corrections (PRs #55 and #57).

# joineRML 0.3.0

## New features

* Added new functions `dynSurv()` and `dynLong()`, which generates survival probabilities and expected longitudinal predictions, respectively, for a new subject conditional on their last measurement time and longitudinal history. Prediction can be implemented using either a first order approximation or a Monte Carlo simulation approach.

* Added an associated `print()` method for `dynSurv` and `dynLong` objects.

* Added an associated `plot()` method for `dynSurv` and `dynLong` objects.

* Added a function `baseHaz()` for extracting the centered and uncentered estimates of the baseline hazard function.

* `print()` and `summary()` now report the total computation time in addition to just the EM algorithm time. This was deemed useful after some examples showed that the time to get initial values was more expensive than the time for the MCEM algorithm to converge.

## Bug patches

* Fixed a bug that prevented models being fitted with no covariates in the survival sub-model, i.e. `Surv() ~ 1`.

* Correction to the vignette description of `mjoint()` arguments.

* Removed `enumintem` package for Sweave vignette to satisfy CRAN checks on macOS (release).

## Housekeeping

* Updated `Makevars` and `Markvars.win` to allow for OpenMP, which can be used by RcppArmadillo.

* Minor tidy-up of in-code comments.

* Minor updates and corrections to documentation.

* Added unit tests for new features.

# joineRML 0.2.2

* This update is an attempt to overcome a FAIL status on the CRAN checks for macOS.

## Housekeeping

* Changed the technical vignette to Rnw (with engine = Sweave) from ltx (with engine = R.rsp) in attempt to remove some CMD check warnings.

* Shortened the vignette to make it compile quicker (removed execution of bootstrapping).

* Lots of tweaks of minor formatting tweaks in the documentation.

# joineRML 0.2.1

## Bug patches

* Fixed a bug that prevented factors with >2 levels being included in the time-to-event sub-model.

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

* Added a check that if multiple repeated longitudinal outcomes were given, that each subject contributes at least one measurement per outcome.

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



