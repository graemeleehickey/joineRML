## Reverse-dependency results

On 24 November 2026, we compared the four current CRAN reverse dependencies using joineRML 0.4.8 and the 0.4.9 candidate on macOS with R 4.5.3. Supporting package versions were identical.

No new problems were found. BCClong 1.0.3 and jmBIG 0.1.4 completed with no errors, warnings, or notes. Broom 1.0.13 had the same spelling NOTE in both runs; its tests had the same unrelated survival deprecation warning and 23 skips for unavailable optional packages. All 39 joineRML-specific expectations passed with both versions, without warnings or skips.

JMbdirect 0.1.1 could not be installed with either version: it imports `survfitjmcs`, which is not exported by FastJM 1.7.2. This prevents a complete compatibility assessment of JMbdirect; the failure is not introduced by joineRML 0.4.9.

Checks included examples, available tests, and R code from vignettes, but did not rebuild PDF manuals or vignette outputs. See [the full report](README.md) for methods, logs, and additional jmBIG prediction checks.
