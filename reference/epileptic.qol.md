# Quality of life data following epilepsy drug treatment

The SANAD (Standard and New Antiepileptic Drugs) study (Marson et al.,
2007) is a randomised control trial of standard and new antiepileptic
drugs, comparing effects on longer term clinical outcomes. Quality of
life (QoL) data were collected by mail at baseline, 3 months, and at 1
and 2 years using validated measures. This data is a subset of the trial
for 544 patients randomised to one of 2 drugs: carbamazepine and
lamotrigine.

## Usage

``` r
data(epileptic.qol)
```

## Format

A data frame with 1853 observations on the following 9 variables:

- `id`:

  patients identifier; in total there are 544 patients.

- `with.time`:

  number of days between registration and the earlier of treatment
  failure or study analysis time.

- `trt`:

  a factor with levels `CBZ` and `LTG` denoting carbamazepine and
  lamotrigine, respectively.

- `with.status`:

  the reason for treatment failure. Coded as `0=`censored;
  `1=`unacceptable adverse effects; `2=`inadequate seizure control.

- `time`:

  the time the quality of life measures were recorded (days). The first
  measurement for each subject is the baseline measurement, however
  there was variability between the time taken to return the
  questionnaires; hence the reason this is non-zero. Similarly, the
  second, third, and fourth follow-up times, which were scheduled for
  3-months, 1-year, and 2-years, respectively, also had variability in
  completion times.

- `anxiety`:

  a continuous measure of anxiety, as defined according to the NEWQOL
  (Newly Diagnosed Epilepsy Quality of Life) assessment. Higher scores
  are indicative of worse QoL.

- `depress`:

  a continuous measure of depression, as defined according to the NEWQOL
  (Newly Diagnosed Epilepsy Quality of Life) assessment. Higher scores
  are indicative of worse QoL.

- `aep`:

  a continuous measure of the Liverpool Adverse Events Profile (AEP), as
  defined according to the NEWQOL (Newly Diagnosed Epilepsy Quality of
  Life) assessment. Higher scores are indicative of worse QoL.

- `with.status2`:

  a binary indicator of composite treatment failure (for any reason),
  coded `status2=1`, or right-censoring `status2=0`.

## Source

SANAD Trial: University of Liverpool. See Jacoby et al. (2015).

## References

Jacoby A, Sudell M, Tudur Smith C, et al. Quality-of-life outcomes of
initiating treatment with standard and newer antiepileptic drugs in
adults with new-onset epilepsy: Findings from the SANAD trial.
*Epilepsia*. 2015; **56(3)**: 460-472.

Marson AG, Appleton R, Baker GA, et al. A randomised controlled trial
examining longer-term outcomes of standard versus new antiepileptic
drugs. The SANAD Trial. *Health Technology Assessment*. 2007;
**11(37)**.

Marson AG, Al-Kharusi AM, Alwaidh M, et al. The SANAD study of
effectiveness of carbamazepine, gabapentin, lamotrigine, oxcarbazepine,
or topiramate for treatment of partial epilepsy: an unblinded randomised
controlled trial. *Lancet*. 2007; **365**: 2007-2013.

Abetz L, Jacoby A, Baker GA, et al. Patient-based assessments of quality
of life in newly diagnosed epilepsy patients: validation of the NEWQOL.
*Epilepsia*. 2000; **41**: 1119-1128.

## See also

[`pbc2`](https://graemeleehickey.github.io/joineRML/reference/pbc2.md),
[`heart.valve`](https://graemeleehickey.github.io/joineRML/reference/heart.valve.md),
[`renal`](https://graemeleehickey.github.io/joineRML/reference/renal.md).
