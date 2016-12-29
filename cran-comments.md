## Resubmission
This is a resubmission owing to an installation ERROR reported by CRAN checks for r-oldrel-windows-ix86+x86_64. In this version I have:

* Patched sigma.mjoint.R to reference S3 method for sigma depending on version number.

## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* local Windows 7 Enterprise SP1 install, R 3.3.1 (via github_install)
* win-builder (devel and release)
* r-oldrel-windows (via R-hub)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

## Other notes

The vignette code can take about 15-20 minutes to run due to the Monte Carlo simulation example.
