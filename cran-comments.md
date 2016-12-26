## Resubmission
This is a resubmission. In this version I have:

* Removed the space from the DOI in the DESCRIPTION.

* Reduced run times for all examples in Rd files. All take approx. 2 seconds on my machine (MacBook Air Intel Core i5, 8GB RAM) and cause no NOTES on the win-builder (devel).

## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* local Windows 7 Enterprise SP1 install, R 3.3.1 (via github_install)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

## Other notes

The vignette code can take about 15-20 minutes to run due to the Monte Carlo simulation example.
