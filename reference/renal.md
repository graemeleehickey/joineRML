# Renal transplantation data

This is a dataset on 407 patients suffering from chronic kidney disease
who underwent a primary renal transplantation with a graft from a
deceased or living donor in the University Hospital of the Catholic
University of Leuven (Belgium) between 21 January 1983 and 16 August
2000. Chronic kidney (renal) disease is a progressive loss of renal
function over a period of months or years through five stages. Each
stage is a progression through an abnormally low and progressively worse
glomerular filtration rate (GFR). The dataset records 3 repeated
measures (2 continuous and 1 binary), and an event time.

## Usage

``` r
data(renal)
```

## Format

This is a list with 4 data frames:

1.  `prot`: repeated measurement data for proteinuria (binary) that
    measures whether the kidneys succeed in sustaining the proteins in
    the blood and not discard them in the urine.

2.  `haem`: repeated measurement data for blood haematocrit level
    (continuous) that measures whether the kidneys produce adequate
    amounts of the hormone erythropoietin that regulates the red blood
    cell production.

3.  `gfr`: repeated measurement data for GFR (continuous) that measures
    the filtration rate of the kidneys.

4.  `surv`: time-to-event data for renal graft failure.

**All datasets** have the common data columns, which are in long format
for the 3 longitudinal data data frames, and 1-per-subject for the
time-to-event data frame:

- `id`:

  number for patient identification.

- `age`:

  age of patient at day of surgery (years).

- `weight`:

  preoperative weight of patient (kg).

- `sex`:

  gender of patient.

- `fuyears`:

  maximum follow up time, with transplant date as the time origin
  (years).

- `failure`:

  censoring indicator (`1=`graft failure and `0=`censored).

**The longitudinal datasets only** contain 2 further columns:

- `time`:

  observed time point, with surgery date as the time origin (years).

- biomarker value:

  a recorded measurement of the biomarker taken at time `time`. The 3
  biomarkers (one per data frame) are:

  - `proteinuria`: recorded as binary indicator: present or not-present.
    Present in the `prot` data.

  - `haematocrit`: recorded as percentage (%) of the ratio of the volume
    of red blood cells to the total volume of blood. Present in the
    `haem` data.

  - `gfr`: measured as ml/min/1.73\\m^2\\. Present in the `gfr` data.

## Source

Dr Dimitris Rizopoulos (<d.rizopoulos@erasmusmc.nl>).

## References

Rizopoulos D, Ghosh, P. A Bayesian semiparametric multivariate joint
model for multiple longitudinal outcomes and a time-to-event. *Stat
Med.* 2011; **30(12)**: 1366-80.

## See also

[`pbc2`](https://graemeleehickey.github.io/joineRML/reference/pbc2.md),
[`heart.valve`](https://graemeleehickey.github.io/joineRML/reference/heart.valve.md),
[`epileptic.qol`](https://graemeleehickey.github.io/joineRML/reference/epileptic.qol.md).
