# joineRML: Joint Modelling of Multivariate Longitudinal Data and Time-to-Event Outcomes

Fits the joint model proposed by Henderson and colleagues (2000)
[doi:10.1093/biostatistics/1.4.465](https://doi.org/10.1093/biostatistics/1.4.465)
, but extended to the case of multiple continuous longitudinal measures.
The time-to-event data is modelled using a Cox proportional hazards
regression model with time-varying covariates. The multiple longitudinal
outcomes are modelled using a multivariate version of the Laird and Ware
linear mixed model. The association is captured by a multivariate latent
Gaussian process. The model is estimated using a Monte Carlo Expectation
Maximization algorithm. This project was funded by the Medical Research
Council (Grant number MR/M013227/1).

## See also

Useful links:

- <https://github.com/graemeleehickey/joineRML>

- Report bugs at <https://github.com/graemeleehickey/joineRML/issues>

## Author

**Maintainer**: Graeme L. Hickey <graemeleehickey@gmail.com>
([ORCID](https://orcid.org/0000-0002-4989-0054))

Authors:

- Pete Philipson <peter.philipson1@newcastle.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-7846-0208))

- Ruwanthi Kolamunnage-Dona <kdrr@liverpool.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-3886-6208))

- Alessandro Gasparini <alessandro.gasparini@ki.se>
  ([ORCID](https://orcid.org/0000-0002-8319-7624))

Other contributors:

- Andrea Jorgensen <aljorgen@liverpool.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-6977-9337)) \[contributor\]

- Paula Williamson <p.r.williamson@liverpool.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-9802-6636)) \[contributor\]

- Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl> (data/renal.rda,
  R/hessian.R, R/vcov.R) \[contributor, data contributor\]

- Medical Research Council (Grant number: MR/M013227/1) \[funder\]
