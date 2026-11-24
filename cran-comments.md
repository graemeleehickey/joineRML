## Other notes

* This is a maintenance release with bug fixes, minor improvements, and an updated minimum R version.

## Test environments

* local macOS (Tahoe 26.5), R 4.5.3
* ubuntu (via GitHub Actions, release + devel)
* macOS (via GitHub Actions, release)
* windows (via GitHub Actions, release)
* windows (via win-builder, old + release + devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

On 24 November 2026, we compared all four current CRAN reverse dependencies using joineRML 0.4.8 and the 0.4.9 candidate on macOS Tahoe 26.7.3 with R 4.5.3 and identical supporting package versions. Checks included examples, available tests, and R code from vignettes; PDF manuals and vignette outputs were not rebuilt.

* No new problems were found. BCClong 1.0.3 and jmBIG 0.1.4 completed with no errors, warnings, or notes.
* Broom 1.0.13 had the same spelling NOTE and unrelated survival deprecation warning in both runs. Its 23 skipped tests concerned unavailable optional packages. All 39 joineRML-specific expectations passed with both versions, without warnings or skips.
* JMbdirect 0.1.1 failed to install with both versions because it imports `survfitjmcs`, which FastJM 1.7.2 does not export. This existing dependency incompatibility prevents a complete assessment of JMbdirect.
