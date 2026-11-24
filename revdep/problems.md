# Problems observed

No new problems were found in the completed checks of BCClong, broom, and jmBIG when comparing joineRML 0.4.8 with 0.4.9. JMbdirect could not be checked beyond installation; see [the failure report](failures.md).

Broom 1.0.13 had the following identical results with both joineRML versions:

- One `R CMD check` NOTE for differences in the spelling-test output, listing words such as `cli`, `namespacing`, and `testthat`.
- Its test suite recorded 3,083 passed expectations, no failures, one warning about the deprecated `survival::Surv(type = "mstate")` argument, and 23 skipped tests for unavailable optional packages.
- Its joineRML-specific tests passed all 39 expectations without warnings or skips in a separate run.

The 22 unavailable optional packages were AUC, binGroup, btergm, drc, ergm, gee, glmnetUtils, gmm, irlba, Kendall, ks, lm.beta, lmodel2, mediation, mfx, modeldata, plm, poLCA, robust, spatialreg, systemfit, and vars. This limits coverage of broom's other model types, but did not prevent its joineRML checks.

BCClong and jmBIG had no errors, warnings, or notes. BCClong's test suite passed 29 expectations without warnings or skips. JMbdirect's installation failure and broom's existing findings are documented in the [CRAN-version logs](logs/cran.log) and [candidate-version logs](logs/candidate.log).
