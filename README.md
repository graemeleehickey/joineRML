
<!-- README.md is generated from README.Rmd. Please edit that file -->

# joineRML <img src="man/figures/hex.png" width = "175" height = "200" align="right" />

[![Travis-CI Build
Status](https://travis-ci.org/graemeleehickey/joineRML.svg?branch=master)](https://travis-ci.org/graemeleehickey/joineRML)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/graemeleehickey/joineRML?branch=master&svg=true)](https://ci.appveyor.com/project/graemeleehickey/joineRML)
<!--[![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)-->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/joineRML)](https://CRAN.R-project.org/package=joineRML)
[![](http://cranlogs.r-pkg.org/badges/joineRML)](https://CRAN.R-project.org/package=joineRML)
[![](https://cranlogs.r-pkg.org/badges/grand-total/joineRML)](https://CRAN.R-project.org/package=joineRML)
[![codecov](https://codecov.io/gh/graemeleehickey/joineRML/branch/master/graph/badge.svg)](https://codecov.io/gh/graemeleehickey/joineRML)
[![Research software
impact](http://depsy.org/api/package/cran/joineRML/badge.svg)](http://depsy.org/package/r/joineRML)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1158231.svg)](https://doi.org/10.5281/zenodo.1158231)

`joineRML` is an extension of the joineR package for fitting joint
models of time-to-event data and multivariate longitudinal data. The
model fitted in joineRML is an extension of the Wulfsohn and Tsiatis
(1997) and Henderson et al. (2000) models, which is comprised of
\((K+1)\)-sub-models: a Cox proportional hazards regression model (Cox,
1972) and a \(K\)-variate linear mixed-effects model - a direct
extension of the Laird and Ware (1982) regression model. The model is
fitted using a Monte Carlo Expectation-Maximization (MCEM) algorithm,
which closely follows the methodology presented by Lin et al. (2002).

## Why use joineRML?

As noted in Hickey et al. (2016), there is a lack of statistical
software available for fitting joint models to multivariate longitudinal
data. This is contrary to a growing methodology in the statistical
literature. `joineRML` is intended to fill this void.

# Example

The main workhorse function is `mjoint`. As a simple example, we use the
`heart.valve` dataset from the package and fit a bivariate joint model.

``` r
library(joineRML)
data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

set.seed(12345)
fit <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    timeVar = "time")
```

The fitted model is assigned to `fit`. We can apply a number of
functions to this object, e.g. `coef`, `logLik`, `plot`, `print`,
`ranef`, `fixef`, `summary`, `AIC`, `getVarCov`, `vcov`, `confint`,
`sigma`, `update`, `formula`, `resid`, and `fitted`. In addition,
several special functions have been added, including `dynSurv`,
`dynLong`, and `baseHaz`, as well as plotting functions for objects
inheriting from the `dynSurv`, `dynLong`, `ranef`, and `mjoint`
functions. For example,

``` r
summary(fit)
plot(fit, param = 'gamma')
```

`mjoint` automatically estimates approximate standard errors using the
empirical information matrix (Lin et al., 2002), but the `bootSE`
function can be used as an alternative.

# Errors and updates

If you spot any errors or wish to see a new feature added, please file
an issue at <https://github.com/graemeleehickey/joineRML/issues> or
email [Graeme Hickey](mailto:graeme.hickey@liverpool.ac.uk).

# Further learning

For an overview of the model estimation being performed, please see the
technical vignette, which can be accessed by

``` r
vignette('technical', package = 'joineRML')
```

For a demonstration of the package, please see the introductory
vignette, which can be accessed by

``` r
vignette('joineRML', package = 'joineRML')
```

# Funding

This project is funded by the [Medical Research
Council](http://www.mrc.ac.uk) (Grant number
MR/M013227/1).

![](http://www.mrc.ac.uk/mrc/includes/themes/MRC/images/template/desktop/logo.png)

# Using the latest developmental version

To install the latest **developmental version**, you will need R version
(version 3.3.0 or higher) and some additional software depending on what
platform you are using.

## Windows

If not already installed, you will need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). Choose the
version that corresponds to the version of R that you are using.

## Mac OSX

If not already installed, you will need to install Xcode Command Line
Tools. To do this, open a new terminal and run

``` bash
$ xcode-select --install
```

## From R

The latest developmental version will not yet be available on CRAN.
Therefore, to install it, you will need `devtools`. You can check you
are using the correct version by running

Once the prerequisite software is installed, you can install `joineRML`
by running the following command in an R console

``` r
library('devtools')
install_github('graemeleehickey/joineRML')
```

# Compatibility with [`broom`](https://github.com/tidymodels/broom/)

Tidiers methods for objects of class `mjoint` (i.e. models fit with
`joineRML`) are included in the
[`broom`](https://github.com/tidymodels/broom/) package; this provides
methods that allow extracting model estimates, predictions, and
comparing models in a straightforward way.

See `vignette(topic = "joineRML-broom", package = "joineRML")` for
further details and examples.

# References

1.  Cox DR. Regression models and life-tables. *J R Stat Soc Ser B Stat
    Methodol.* 1972; **34(2)**: 187-220.

2.  Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
    measurements and event time data. *Biostatistics.* 2000; **1(4)**:
    465-480.

3.  Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. Joint
    modelling of time-to-event and multivariate longitudinal outcomes:
    recent developments and issues. *BMC Med Res Methodol.* 2016;
    **16(1)**: 117.

4.  Laird NM, Ware JH. Random-effects models for longitudinal data.
    *Biometrics.* 1982; **38(4)**: 963-974.

5.  Lin H, McCulloch CE, Mayne ST. Maximum likelihood estimation in the
    joint analysis of time-to-event and multiple longitudinal variables.
    *Stat Med.* 2002; **21**: 2369-2382.

6.  Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
    data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.
