## Other notes

* Addresses all NOTES from CRAN Package Check Results.
* I have stopped several examples from running in order to avoid the NOTE on the Debian system. I cannot reproduce these, but am happy to wrap in `\dontrun{}` blocks.

## Test environments

* local macOS (Sonoma 14.6.1) install, R 4.4.2
* ubuntu (via GitHub actions, release + devel)
* macOS (via GitHub actions, release)
* windows (via appveyor CI, release)
* windows (via GitHub actions, release)
* windows (via win-builder, old + release + devel)

## R CMD check results

0 errors | 0 warnings | 2 notes

Win-Builder NOTE: "checking CRAN incoming feasibility ... NOTE"

Debian NOTE: "Re-building vignettes had CPU time 5.5 times elapsed time" -- I cannot reproduce this. It passes every other platform. I suspect the qpdf compression tool is the cause, but that was requested by CMD tests.

## Reverse dependencies

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
