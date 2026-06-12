# The baseline hazard estimate of an `mjoint` object

This function returns the (baseline) hazard increment from a fitted
`mjoint` object. In addition, it can report either the *uncentered* or
the more ubiquitous *centered* version.

## Usage

``` r
baseHaz(object, centered = TRUE, se = FALSE)
```

## Arguments

- object:

  an object inheriting from class `mjoint` for a joint model of
  time-to-event and multivariate longitudinal data.

- centered:

  logical: should the baseline hazard be for the mean-centered
  covariates model or not? Default is `centered=TRUE`. See **Details**.

- se:

  logical: should standard errors be approximated for the hazard
  increments? Default is `se=FALSE`.

## Value

A `data.frame` with two columns: the unique failure times and the
estimate baseline hazard. If `se=TRUE`, then a third column is appended
with the corresponding standard errors (for the centred case).

## Details

When covariates are included in the time-to-event sub-model,
[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
automatically centers them about their respective means. This also
applies to non-continuous covariates, which are first coded using a
dummy-transformation for the design matrix and subsequently centered.
The reason for the mean-centering is to improve numerical stability, as
the survival function involves exponential terms. Extracting the
baseline hazard increments from
[`mjoint.object`](https://graemeleehickey.github.io/joineRML/reference/mjoint.object.md)
returns the Breslow hazard estimate (Lin, 2007) that corresponds to this
mean-centered model. This is the same as is done in the R `survival`
package when using
[`coxph.detail`](https://rdrr.io/pkg/survival/man/coxph.detail.html)
(Therneau and Grambsch, 2000). If the user wants to access the baseline
hazard estimate for the model in which no mean-centering is applied,
then they can use this function, which scales the mean-centered baseline
hazard by

\$\$\exp\\-\bar{w}^\top \gamma_v\\,\$\$

where \\\bar{w}\\ is a vector of the means from the time-to-event
sub-model design matrix.

## References

Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox
Model.* New Jersey: Springer-Verlag; 2000.

Lin DY. On the Breslow estimator. *Lifetime Data Anal.* 2007; **13(4)**:
471-480.

## See also

[`mjoint`](https://graemeleehickey.github.io/joineRML/reference/mjoint.md)
and [`coef`](https://rdrr.io/r/stats/coef.html).

## Author

Graeme L. Hickey (<graemeleehickey@gmail.com>)

## Examples

``` r

if (FALSE) { # \dontrun{
# Fit a joint model with bivariate longitudinal outcomes

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$log.grad) & !is.na(heart.valve$log.lvmi), ]

fit2 <- mjoint(
    formLongFixed = list("grad" = log.grad ~ time + sex + hs,
                         "lvmi" = log.lvmi ~ time + sex),
    formLongRandom = list("grad" = ~ 1 | num,
                          "lvmi" = ~ time | num),
    formSurv = Surv(fuyrs, status) ~ age,
    data = list(hvd, hvd),
    inits = list("gamma" = c(0.11, 1.51, 0.80)),
    timeVar = "time",
    verbose = TRUE)
baseHaz(fit2, centered = FALSE)
} # }
```
