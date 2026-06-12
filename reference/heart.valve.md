# Aortic valve replacement surgery data

This is longitudinal data on an observational study on detecting effects
of different heart valves, differing on type of tissue, implanted in the
aortic position. The data consists of longitudinal measurements on three
different heart function outcomes, after surgery occurred. There are
several baseline covariates available, and also survival data.

## Usage

``` r
data(heart.valve)
```

## Format

This is a data frame in the unbalanced format, that is, with one row per
observation. The data consists in columns for patient identification,
time of measurements, longitudinal multiple longitudinal measurements,
baseline covariates, and survival data. The column names are identified
as follows:

- `num`:

  number for patient identification.

- `sex`:

  gender of patient (`0=`Male and `1=`Female).

- `age`:

  age of patient at day of surgery (years).

- `time`:

  observed time point, with surgery date as the time origin (years).

- `fuyrs`:

  maximum follow up time, with surgery date as the time origin (years).

- `status`:

  censoring indicator (`1=`died and `0=`lost at follow up).

- `grad`:

  valve gradient at follow-up visit.

- `log.grad`:

  natural log transformation of `grad`.

- `lvmi`:

  left ventricular mass index (standardised) at follow-up visit.

- `log.lvmi`:

  natural log transformation of `lvmi`.

- `ef`:

  ejection fraction at follow-up visit.

- `bsa`:

  preoperative body surface area.

- `lvh`:

  preoperative left ventricular hypertrophy.

- `prenyha`:

  preoperative New York Heart Association (NYHA) classification
  (`1=`I/II and `3=`III/IV).

- `redo`:

  previous cardiac surgery.

- `size`:

  size of the valve (millimeters).

- `con.cabg`:

  concomitant coronary artery bypass graft.

- `creat`:

  preoperative serum creatinine (\\\mu\\mol/mL).

- `dm`:

  preoperative diabetes.

- `acei`:

  preoperative use of ace inhibitor.

- `lv`:

  preoperative left ventricular ejection fraction (LVEF) (`1=`good,
  `2=`moderate, and `3=`poor).

- `emergenc`:

  operative urgency (`0=`elective, `1 = `urgent, and `3=`emergency).

- `hc`:

  preoperative high cholesterol (`0=`absent, `1 =`present treated, and
  `2=`present untreated).

- `sten.reg.mix`:

  aortic valve haemodynamics (`1=`stenosis, `2=`regurgitation,
  `3=`mixed).

- `hs`:

  implanted aortic prosthesis type (`1=`homograft and `0=`stentless
  porcine tissue).

## References

Lim E, Ali A, Theodorou P, Sousa I, Ashrafian H, Chamageorgakis T,
Duncan M, Diggle P, Pepper J. A longitudinal study of the profile and
predictors of left ventricular mass regression after stentless aortic
valve replacement. *Ann Thorac Surg.* 2008; **85(6)**: 2026-2029.

## See also

[`pbc2`](https://graemeleehickey.github.io/joineRML/reference/pbc2.md),
[`renal`](https://graemeleehickey.github.io/joineRML/reference/renal.md),
[`epileptic.qol`](https://graemeleehickey.github.io/joineRML/reference/epileptic.qol.md).
