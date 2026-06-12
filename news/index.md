# Changelog

## joineRML 0.4.7.9000

### Bug fixes

- [`bootSE()`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
  parallel branch condition was `ncores >= 1L` (always `TRUE`), making
  the serial code path with the progress bar unreachable. Corrected to
  `ncores > 1L`.

- [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  now throttles the number of threads used by `RcppArmadillo` during the
  MCEM algorithm (restored on exit), so the previously imported
  `armadillo_throttle_cores()` and `armadillo_reset_cores()` functions
  are actually invoked.

- [`simData()`](https://graemeleehickey.github.io/joineRML/reference/simData.md)
  now uses `drop = FALSE` when subsetting the random effects matrices in
  the Gompertz (`model = "intslope"`) survival-time calculation,
  preventing a latent dimension-dropping error when only a single
  subject remains in the branch.

- Internal initial-value calculation for the survival sub-model
  (`initsSurv()`) no longer relies on a leaked loop variable when
  indexing `timeVar`; it now uses the (outcome-invariant) first time
  variable explicitly in the balanced-data case.

### Housekeeping

- Added test coverage for the dynamic prediction functions
  ([`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
  and
  [`dynLong()`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md),
  including first-order and simulated prediction, `u`/`horizon`
  arguments, input validation, and print methods), their plot methods
  ([`plot.dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/plot.dynSurv.md)
  and
  [`plot.dynLong()`](https://graemeleehickey.github.io/joineRML/reference/plot.dynLong.md)),
  `print.summary.mjoint()`,
  [`baseHaz()`](https://graemeleehickey.github.io/joineRML/reference/baseHaz.md),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and several
  accessor methods.

- Added a pkgdown site configuration, GitHub Actions workflow for
  publishing the site to GitHub Pages, and the pkgdown site URL to
  `DESCRIPTION`. The Articles menu links to the technical details
  vignette PDF.

- Added `pkgdown` to `Suggests`.

- Updated GitHub Actions workflows to use `actions/checkout@v5`
  (Node.js 24) and modernised the test-coverage workflow to the current
  r-lib template.

- Added alt text to README images and download badges so the pkgdown
  home page passes accessibility checks while retaining the MRC logo.

- Removed the AppVeyor CI badge and references (`appveyor.yml`
  `.Rbuildignore` entry, `skip_on_appveyor()` in tests, and the
  cran-comments test environment), as the project now relies on GitHub
  Actions.

- Bumped the minimum R version to 4.1.0.

- Added `.posit` to `.Rbuildignore` so local Positron metadata is not
  included in package builds.

- Migrated the test suite to the testthat 3rd edition (added
  `Config/testthat/edition: 3` and `testthat (>= 3.0.0)`), removing
  deprecated `context()` calls, replacing `expect_is()` and
  [`is.ggplot()`](https://ggplot2.tidyverse.org/reference/is_tests.html)
  with `expect_s3_class()`, `expect_type()`, and base class checks, and
  narrowing the unmatched-`inits` warning test to assert the intended
  warning only.

- [`plot.ranef.mjoint()`](https://graemeleehickey.github.io/joineRML/reference/plot.ranef.mjoint.md)
  now uses `geom_errorbar(orientation = "y")` instead of deprecated
  [`geom_errorbarh()`](https://ggplot2.tidyverse.org/reference/geom_linerange.html),
  eliminating ggplot2 4.0 deprecation warnings.

- Removed obsolete `bindrcpp` from `Suggests` in `DESCRIPTION`; the
  package has been superseded and is not used anywhere in the codebase.

- Replaced `1:n`-style sequences with
  [`seq_len()`](https://rdrr.io/r/base/seq.html) and
  [`seq_along()`](https://rdrr.io/r/base/seq.html) in iteration paths
  that could otherwise iterate incorrectly over a zero-length object.

- Bumped `RoxygenNote` to 7.3.3 to match the installed version.

- [`simData()`](https://graemeleehickey.github.io/joineRML/reference/simData.md)
  now reports the simulated event rate via
  [`message()`](https://rdrr.io/r/base/message.html) rather than
  [`cat()`](https://rdrr.io/r/base/cat.html), so the output can be
  suppressed.

- [`summary()`](https://rdrr.io/r/base/summary.html) documentation in
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  now refers to the `bootSE` argument by name
  (`summary(fit_obj, bootSE = boot_obj)`).

- [`vcov.mjoint()`](https://graemeleehickey.github.io/joineRML/reference/vcov.mjoint.md)
  now warns when it falls back to the Moore-Penrose pseudo-inverse, as
  the approximate standard errors may be unreliable in that case.

- Fixed documentation typos: “left hand-hand side” in
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md),
  “comprised on” in the package description, and “The choice os” in
  [`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
  and
  [`dynLong()`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md).

- Removed empty test body for `"argument not a summary.mjoint object"`
  in `test-inputs.R`.

## joineRML 0.4.7

CRAN release: 2025-02-04

### Housekeeping

- Updated README: badges and MRC logo.

- Removed CXX11 requirement from Makevars files.

- Update
  [`plot.ranef.mjoint()`](https://graemeleehickey.github.io/joineRML/reference/plot.ranef.mjoint.md)
  function due to deprecated
  [`aes_string()`](https://ggplot2.tidyverse.org/reference/aes_.html)
  call.

- Fixed deprecated documentation for `joineRML.R`.

- Updated broken links in `plot.ranef.mjoint.R` and `vcov.mjoint.R`.

- Updated GitHub actions workflows.

- Ran reverse dependency checks using `revdepcheck` package.

- Updated CITATION file due to deprecated function and typo.

- Wrapped a bunch of examples from being run on CRAN submission to avoid
  Debian system NOTE.

- Reduced computational time for vignettes to to try and avoid Debian
  system NOTE on CRAN submission.

- Add R-hub workflow to enable checks for platforms not covered already.

- Imported `RcppArmadillo` functions `armadillo_throttle_cores` and
  `armadillo_reset_cores` to try and prevent CRAN R CMD checks using
  multiple threads.

## joineRML 0.4.6

CRAN release: 2023-01-20

### Housekeeping

- Package SUGGESTS are now conditional in the vignettes per CRAN
  requirement.

- [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  now returns the input `inits` and the interim calculated initial
  values calculated for the EM algorithm step: `inits.long` and
  `inits.surv`. Requested by James Murray.

- Updated Makevars and Makevars.win due to compilation issues.

## joineRML 0.4.5

CRAN release: 2021-01-05

### Housekeeping

- Updated Makevars and Makevars.win to remove clang4 dependency.

- Updated URL in DESCRIPTION.

- Fixed couple of typos.

- Moved to GitHub Actions CI.

- Updated URLs.

## joineRML 0.4.4

CRAN release: 2020-04-09

### New features

- Added `broom` compatibility.

### Bug patches

- Fixed issue that required subject IDs to be first column in dataset.

### Housekeeping

- Added [@ellessenne](https://github.com/ellessenne) as package author.

- Updates to DESCRIPTION.

- Documentation updates.

## joineRML 0.4.3

CRAN release: 2020-02-17

### Housekeeping

- Graeme Hickey taken over the package as creator and maintainer.

- Updates to DESCRIPTION

- Updates to Makevars and Makevars.win

- Updates to appveyor.yml and .travis.yml settings

- Updates to .gitignore

## joineRML 0.4.2

CRAN release: 2018-05-28

### Bug patches

- Fixed issue with citation date due to DATE not being present in the
  DESCRIPTION file.

### Housekeeping

- Added Zenodo DOI badge to README.

- Added ORCID IDs for authors.

- Remove unused objects once finished with them to try and reduce memory
  overheads.

- Changed maintainer to Pete Philipson.

## joineRML 0.4.1

CRAN release: 2018-01-21

### New features

- Added smoothed predicted survival curves to the `plot.dynSuv()`.
  Smoothing is based on the constrained B-splines method.

- [`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
  now includes an argument to specify a horizon time from the last known
  observation time.

- [`simData()`](https://graemeleehickey.github.io/joineRML/reference/simData.md)
  now includes an argument to choose multivariate *t*-distributed random
  effects with varying degrees of freedom, thus allowing for sensitivity
  analyses of heavier tail distributions.

### Housekeeping

- Minor corrections to documentation.

- Minor bug fixes to plotting functions.

- Added MRC to DESCRIPTION as funder.

- Added Depsy badge to README.

## joineRML 0.4.0

CRAN release: 2017-11-11

### New features

- The MCEM algorithm can be implemented with vanilla Monte Carlo and
  quasi-Monte Carlo (using either the scrambled Sobol sequence or the
  Halton sequence). This is implemented through the `type` argument
  passed to the list of `control` parameters in
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md).

- [`bootSE()`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
  now has option to use parallel computing via the `foreach` package.

### Bug patches

- Fixed some small errors in the `epileptic.qol` dataset.

- Fixed situation where a tibble might be given as the dataset
  ([\#55](https://github.com/graemeleehickey/joineRML/issues/55)
  [@ellessenne](https://github.com/ellessenne)).

- Catches errors in bootstrap due to “bad” data and automatically
  restarts the bootstrap
  ([\#57](https://github.com/graemeleehickey/joineRML/issues/57)
  [@ellessenne](https://github.com/ellessenne))

### Housekeeping

- Added hex sticker badge.

- Moved the make files and raw data for `qol` and `renal` datasets into
  the `~/data-raw/` directory.

- Added [@ellessenne](https://github.com/ellessenne) as package
  contributor for testing and bug corrections (PRs
  [\#55](https://github.com/graemeleehickey/joineRML/issues/55) and
  [\#57](https://github.com/graemeleehickey/joineRML/issues/57)).

## joineRML 0.3.0

CRAN release: 2017-07-23

### New features

- Added new functions
  [`dynSurv()`](https://graemeleehickey.github.io/joineRML/reference/dynSurv.md)
  and
  [`dynLong()`](https://graemeleehickey.github.io/joineRML/reference/dynLong.md),
  which generates survival probabilities and expected longitudinal
  predictions, respectively, for a new subject conditional on their last
  measurement time and longitudinal history. Prediction can be
  implemented using either a first order approximation or a Monte Carlo
  simulation approach.

- Added an associated [`print()`](https://rdrr.io/r/base/print.html)
  method for `dynSurv` and `dynLong` objects.

- Added an associated
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
  `dynSurv` and `dynLong` objects.

- Added a function
  [`baseHaz()`](https://graemeleehickey.github.io/joineRML/reference/baseHaz.md)
  for extracting the centered and uncentered estimates of the baseline
  hazard function.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) now report the
  total computation time in addition to just the EM algorithm time. This
  was deemed useful after some examples showed that the time to get
  initial values was more expensive than the time for the MCEM algorithm
  to converge.

### Bug patches

- Fixed a bug that prevented models being fitted with no covariates in
  the survival sub-model, i.e. `Surv() ~ 1`.

- Correction to the vignette description of
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  arguments.

- Removed `enumintem` package for Sweave vignette to satisfy CRAN checks
  on macOS (release).

### Housekeeping

- Updated `Makevars` and `Markvars.win` to allow for OpenMP, which can
  be used by RcppArmadillo.

- Minor tidy-up of in-code comments.

- Minor updates and corrections to documentation.

- Added unit tests for new features.

## joineRML 0.2.2

CRAN release: 2017-05-01

- This update is an attempt to overcome a FAIL status on the CRAN checks
  for macOS.

### Housekeeping

- Changed the technical vignette to Rnw (with engine = Sweave) from ltx
  (with engine = R.rsp) in attempt to remove some CMD check warnings.

- Shortened the vignette to make it compile quicker (removed execution
  of bootstrapping).

- Lots of tweaks of minor formatting tweaks in the documentation.

## joineRML 0.2.1

CRAN release: 2017-04-25

### Bug patches

- Fixed a bug that prevented factors with \>2 levels being included in
  the time-to-event sub-model.

- Fixed package registration, which strangely broke on R 3.4.0 for OSX
  platform.

### Housekeeping

- `joineRML` version 0.2.1 will depend on R version \>=3.3.0 to remedy
  issue with
  [`sigma.mjoint()`](https://graemeleehickey.github.io/joineRML/reference/sigma.mjoint.md)
  S3 issue.

- Added a new badge to the README.

## joineRML 0.2.0

CRAN release: 2017-03-27

### New features

- Add `residual()` and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) functions for
  `mjoint` objects.

- Added a plot function –
  [`plot.ranef.mjoint()`](https://graemeleehickey.github.io/joineRML/reference/plot.ranef.mjoint.md)
  – for `ranef.mjoint` objects.

- [`bootSE()`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
  now uses control parameters from fitted model and allows for
  individual parameter overwrite.

- Added a check that any initial covariance matrix given is
  positive-definite.

- Added a check that dimensions of any `inits` given match the
  calculated dimensions from the model formulae.

- Added a check that if multiple repeated longitudinal outcomes were
  given, that each subject contributes at least one measurement per
  outcome.

- Changed the API so that `postRE` and `approx.se` arguments are now
  replaced by a single `pfs` argument, which stands for post fit
  statistics. If `TRUE`, then *both* the approximate SEs and the BLUPs
  (and SEs) are calculated and returned. This change is to facilitate
  other post fit statistics, e.g. residuals.

- [`sampleData()`](https://graemeleehickey.github.io/joineRML/reference/sampleData.md)
  now allows for sampling *without* replacement.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`plotConvergence()`](https://graemeleehickey.github.io/joineRML/reference/plotConvergence.md)
  now have the option to discard burn-in phase iterations from the MCEM
  algorithm.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`plotConvergence()`](https://graemeleehickey.github.io/joineRML/reference/plotConvergence.md)
  now plot the log-likelihood trace.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) now calculates the
  variance-covariance matrix rather than inside
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  and then extracting it. It also utilises the QR-decomposition inverse
  function and the Moore-Penrose matrix inverse, as in some cases the
  matrix was nearly singular.

- `hessian()` (and therefore
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)) now calculate the
  contribution for the random effect variance terms rather than the
  random effect precision (1 divided by the variance) terms. The correct
  contribution to the score for off-diagonal terms is now also impleted.

- The left-hand side of the `formLongFixed` now handles transformations.

### Bug patches

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) now returns the
  variance-covariance matrix as intended. Previously it was only
  returning the profile empirical information matrix.

- Patched a major bug in `gammaUpdate()` where ties in failure times
  were not being properly handled. The code for `gammaUpdate_approx()`
  always worked fine, as it was based only on the score vector. This bug
  manifested when
  [`bootSE()`](https://graemeleehickey.github.io/joineRML/reference/bootSE.md)
  was called due to the resampling with replacement yielding datasets
  with many more ties than in the original dataset used to fit the
  model. To fix it, the information matrix required scaling at each
  failure time by the number of failures in the data. The formula for
  the information matrix in the Technical Details vignette has also been
  updated.

- Patched a minor bug with prevented univariate random-intercept models
  for being fitted.

- Patched a minor bug in plotting convergence traces.

- Patched a minor bug with bootstrapping univariate joint models without
  passing the MLEs as the initial values to the
  [`mjoint()`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
  call.

### Housekeeping

- Renamed `approxSE()` function to `hessian()`.

- Renamed `control` argument `earlyPhase` to `burnin`.

- Default settings for `control` parameters updated based on practical
  experience.

- Package now `Depends` on `survival` and `nlme` rather than `Imports`
  to allow [`require()`](https://rdrr.io/r/base/library.html) statements
  to be removed from code.

- Prevented `roxygen` from exporting all functions.

- Fixed Imports following CRAN Checks of v0.1.1.

- Minor documentation edits and corrections.

- Minor code tidy-up with slight speed-ups and stabilisations where
  found.

- LICENSE upgraded to GPL-3 to be compatible with `joineR` v1.1.0.

- Removed internal function `fast_nearPD()` from package code as was
  unused.

- Removed internal function `EexpArma()` from package code as was
  unused.

- Added unit tests.

- Registered native C++ routines and disabled symbol search to satisfy
  CRAN CMD checks.

## joineRML 0.1.1

CRAN release: 2016-12-30

- Patched sigma.R roxygen documentation to handle sigma S3 method
  changing from lme4 to stats in R v3.3.0.

## joineRML 0.1.0

CRAN release: 2016-12-27

- First package release.
