# Reverse-dependency checks

Graeme L. Hickey · 24 November 2026

No new problems were found when comparing CRAN joineRML 0.4.8 with the 0.4.9 candidate at commit `de53539b81355654b69ce2d0eda2d3b9cc0183c0`. Three reverse dependencies completed their checks. JMbdirect could not be installed with either version because of an existing incompatibility with FastJM, so its compatibility with joineRML 0.4.9 remains unverified.

## Results

The four direct reverse dependencies were identified from CRAN's source package index, including `Depends`, `Imports`, `LinkingTo`, and `Suggests`.

| Package | Version | Uses joineRML through | CRAN 0.4.8 | Candidate 0.4.9 |
|:--------|:--------|:---------------------|:-----------|:----------------|
| BCClong | 1.0.3 | Suggests | OK | OK |
| broom | 1.0.13 | Suggests | 1 NOTE | 1 NOTE |
| jmBIG | 0.1.4 | Imports | OK | OK |
| JMbdirect | 0.1.1 | Imports | Installation ERROR | Installation ERROR |

These are the results recorded in `00check.log`, including broom's spelling NOTE. The summaries returned by `rcmdcheck` omitted that NOTE. See [existing problems](problems.md) and the [installation failure](failures.md) for details.

## Environment and scope

| Setting | Value |
|:--------|:------|
| R | 4.5.3 (2026-03-11) |
| Platform | macOS Tahoe 26.7.3, aarch64-apple-darwin20 |
| Locale | C.UTF-8 |
| Time zone | Europe/London |
| Check command | `R CMD check --no-manual --no-build-vignettes` |
| Check runner | rcmdcheck 1.4.0 |
| Test settings | `NOT_CRAN=true`, `_R_CHECK_FORCE_SUGGESTS_=false` |

The two joineRML versions were installed from source into separate temporary libraries. Supporting packages were shared between runs; comparison of the effective libraries confirmed that joineRML was the only version difference among 992 available packages. Their versions are recorded in [dependency-versions.csv](dependency-versions.csv).

Examples, available tests, and R code from vignettes were checked. PDF manuals and rendered vignette outputs were not rebuilt. Broom had 22 unavailable optional packages and skipped 23 tests for those packages in each run; all its joineRML-specific tests ran. These local checks do not establish results on Windows or Linux.

## Additional checks of joineRML use

- Broom's `test-joineRML.R` passed all 39 expectations with each version: no failures, warnings, errors, or skips. This covers its `tidy()`, `glance()`, and `augment()` checks, including bootstrap standard errors and expected errors.
- `jmBIG::predJRML()` completed with both versions for the same saved univariate and bivariate heart-valve fits. All longitudinal summaries were finite, and survival summaries were finite and between 0 and 1. The calls used subject 10, measurements through year 5, prediction at year 8, `M = 50`, and seed 20261124. Numerical equality was not required because the prediction sampler changed in 0.4.9. These checks assess whether the functions work together, not prediction accuracy or interval coverage.

The fitted examples used complete log-gradient and log-LVMI records with `num <= 50` from `heart.valve` (43 subjects, eight events, and 144 observations per outcome). Both used `Surv(fuyrs, status) ~ age`. The univariate model used `log.grad ~ time + sex + hs` with a random intercept; the bivariate model added `log.lvmi ~ time + sex` with a random intercept and slope. Fitting used `timeVar = "time"`, `type = "halton"`, `nMC = 256`, `nMCmax = 1024`, `burnin = 10`, `mcmaxIter = 100`, `convCrit = "abs"`, `tol0 = 0.02`, and seed 20261124. Initial `gamma` values were `c(0.11, 1.51)` and `c(0.11, 1.51, 0.80)`, respectively. The same fits were reused in both comparisons.

## Check settings and evidence

[CRAN-version logs](logs/cran.log) and [candidate-version logs](logs/candidate.log) retain all package check logs, the JMbdirect installation errors, the test output for BCClong and broom, and the additional checks. Temporary workspace paths are replaced by `<WORKDIR>`; the initial file-copy messages from broom's focused checks are omitted. [Source checksums](sources.csv) identify the package archives used.

The package checks used the following settings in separate R sessions, with the appropriate joineRML library first in `.libPaths()`, followed by the shared supporting libraries. Source archives were the versions listed above.

```r
check_env <- c(
  NOT_CRAN = "true",
  `_R_CHECK_FORCE_SUGGESTS_` = "false",
  `_R_CHECK_CRAN_INCOMING_REMOTE_` = "false",
  RGL_USE_NULL = "true",
  OMP_THREAD_LIMIT = "1",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  R_INSTALL_NCPUS = "1",
  TZ = "Europe/London"
)
rcmdcheck::rcmdcheck(
  path = source_archive,
  args = c("--no-manual", "--no-build-vignettes"),
  check_dir = output_directory,
  libpath = .libPaths(),
  quiet = TRUE,
  timeout = 1800,
  error_on = "never",
  env = check_env
)
```

The focused broom checks used a separate copy of its source `tests/testthat` directory, including helpers and fixtures:

```r
Sys.setenv(NOT_CRAN = "true")
testthat::test_dir(
  path = broom_test_directory,
  filter = "joineRML",
  package = "broom",
  load_package = "installed",
  reporter = "summary",
  stop_on_failure = FALSE
)
```

For each saved fit, the jmBIG calls used:

```r
nd <- subset(fit$data[[1]], num == 10 & time <= 5)
nd$id <- nd$num
sd <- fit$survData
sd$id <- sd$num
set.seed(20261124)
pred <- jmBIG::predJRML(
  model = list(pseudoMod = fit), ids = 10,
  dtlong = nd, dtsurv = sd, u = 8, M = 50, progress = FALSE
)
stopifnot(
  all(vapply(pred$plong[[1]]$pred,
             function(x) all(is.finite(as.matrix(x))), logical(1))),
  all(is.finite(as.matrix(pred$psurv[[1]]$pred))),
  all(as.matrix(pred$psurv[[1]]$pred[, -1]) >= 0),
  all(as.matrix(pred$psurv[[1]]$pred[, -1]) <= 1)
)
```
