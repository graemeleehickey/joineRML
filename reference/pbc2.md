# Mayo Clinic primary biliary cirrhosis data

This data is from the Mayo Clinic trial in primary biliary cirrhosis
(PBC) of the liver conducted between 1974 and 1984. A total of 424 PBC
patients, referred to Mayo Clinic during that ten-year interval met
eligibility criteria for the randomized placebo controlled trial of the
drug D-penicillamine, but only the first 312 cases in the data set
participated in the randomized trial. Therefore, the data here are for
the 312 patients with largely complete data.

## Usage

``` r
data(pbc2)
```

## Format

A data frame with 1945 observations on the following 20 variables:

- `id`:

  patients identifier; in total there are 312 patients.

- `years`:

  number of years between registration and the earlier of death,
  transplantation, or study analysis time.

- `status`:

  a factor with levels `alive`, `transplanted` and `dead`.

- `drug`:

  a factor with levels `placebo` and `D-penicil`.

- `age`:

  at registration in years.

- `sex`:

  a factor with levels `male` and `female`.

- `year`:

  number of years between enrollment and this visit date, remaining
  values on the line of data refer to this visit.

- `ascites`:

  a factor with levels `No` and `Yes`.

- `hepatomegaly`:

  a factor with levels `No` and `Yes`.

- `spiders`:

  a factor with levels `No` and `Yes`.

- `edema`:

  a factor with levels `No edema` (i.e. no edema and no diuretic therapy
  for edema), `edema no diuretics` (i.e. edema present without
  diuretics, or edema resolved by diuretics), and
  `edema despite diuretics` (i.e. edema despite diuretic therapy).

- `serBilir`:

  serum bilirubin in mg/dl.

- `serChol`:

  serum cholesterol in mg/dl.

- `albumin`:

  albumin in mg/dl.

- `alkaline`:

  alkaline phosphatase in U/liter.

- `SGOT`:

  SGOT in U/ml.

- `platelets`:

  platelets per cubic ml/1000.

- `prothrombin`:

  prothrombin time in seconds.

- `histologic`:

  histologic stage of disease.

- `status2`:

  a numeric vector with the value 1 denoting if the patient was dead,
  and 0 if the patient was alive or transplanted.

## Source

[`pbc2`](https://rdrr.io/pkg/JM/man/pbc.html) and
[`pbc`](https://rdrr.io/pkg/survival/man/pbc.html).

## References

Fleming T, Harrington D. *Counting Processes and Survival Analysis*.
1991; New York: Wiley.

Therneau T, Grambsch P. *Modeling Survival Data: Extending the Cox
Model*. 2000; New York: Springer-Verlag.

## See also

[`heart.valve`](https://graemeleehickey.github.io/joineRML/reference/heart.valve.md),
[`renal`](https://graemeleehickey.github.io/joineRML/reference/renal.md),
[`epileptic.qol`](https://graemeleehickey.github.io/joineRML/reference/epileptic.qol.md).
