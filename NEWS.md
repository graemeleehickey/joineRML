# joineRML 0.1.1.9000 (devel version).

## New features

* `bootSE()` now uses control parameters from fitted model.

* Added a check that any initial covariance matrix given is positive-definite.

## Bug patches

* Patched a small bug with prevented univariate random-intercept models for being fitted.

* Patched a small bug in plotting convergence traces.

## Housekeeping

* Package now `Depends` on `survival` and `nlme` rather than `Imports` to allow `require()` statements to be removed from code.

* Prevented `roxygen` from exporting all functions. 

* Fixed Imports following CRAN Checks of v0.1.1.

* Minor documentation edits.

* LICENSE upgraded to GPL-3 to be compatible with `joineR` v1.1.0.

* Removed `fast_nearPD()` from package code as was unused.

* Added unit tests.

# joineRML 0.1.1

* Patched sigma.R roxygen documentation to handle sigma S3 method changing from lme4 to stats in R v3.3.0.

# joineRML 0.1.0

* First package release.



