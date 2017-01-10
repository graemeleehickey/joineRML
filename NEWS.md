# joineRML 0.1.1.9000

* Developmental version of package.

* Added unit tests.

* Patched a small bug with prevented univariate random-intercept models for being fitted.

* Added a check that any initial covariance matrix given is positive-definite.

* LICENSE upgraded to GPL-3 to be compatible with `joineR` v1.1.0.

* Package now `Depends` on `survival` and `nlme` rather than `Imports` to allow `require()` statements to be removed from code.

* Prevented `roxygen` from exporting all functions. 

* Fixed Imports following CRAN Checks of v0.1.1.

* Minor documentation edits.

# joineRML 0.1.1

* Patched sigma.R roxygen documentation to handle sigma S3 method changing from lme4 to stats in R v3.3.0.

# joineRML 0.1.0

* First package release.



