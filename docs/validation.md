# Validation and interpretation

geoCensor separates software checks from validation of statistical inference.

## Software verification

The published 1.2.0 source package passed R CMD check and 1,298 test expectations.
The released commit was checked on macOS, Windows, Linux release, Linux oldrel
and Linux R-devel. The [release workflow](https://github.com/diskoisland/geoCensor/actions/runs/41904919174)
and downloadable [validation receipt](https://github.com/diskoisland/geoCensor/releases/download/v1.2.0/geoCensor_1.2.0_validation.json)
record the exact release source and platform checks. Current development checks
are available from the [workflow page](https://github.com/diskoisland/geoCensor/actions/workflows/R-CMD-check.yaml).

Checks cover parsing and limit recovery, source-text preservation, joins,
unavailable estimates, declustering, method/support audits and installed workflows.
The sediment showcase uses the released analysis functions; the figures are
additional repository documentation. Native-input checks compare selection,
source-record indices and results with separate reference calculations. Software
success does not establish representativeness, measurement accuracy or calibrated
spatial interval coverage.

## Statistical boundaries

| Output | Interpretation |
| --- | --- |
| Censored lognormal geometric mean | Working-model description of the supplied usable records; model adequacy and target population need review. |
| Ordinary model interval | Requires the distributional model and independent-observation assumptions; spatial dependence is not accounted for. |
| Declustering sensitivity | Changes the descriptive influence of occupied cells at stated scales; does not reconstruct unobserved terrain or fully correct preferential sampling. |
| Weight effective sample size | Describes weight imbalance; it is not a count of independent locations. |
| Empirical CDF bounds | Possible cumulative fractions conditional on reported values/limits and weights; these are not confidence bands. |
| Experimental spatial intervals | Restricted to research and explicit validation; not promoted as calibrated for routine NGDB inference. |
| All-left or otherwise unavailable fit | Retained as unavailable; no concentration is imputed. |

Routine declustering-weighted confidence limits remain absent. Simulation work
has exposed undercoverage in tested spatial settings, including controls using
covariance parameters estimated from the observations. Further covariance-refitting,
censoring and model-misspecification experiments are required before promoting
spatial inference. Numerical convergence and a high number of bootstrap fits do
not establish coverage.

## Interpreting environmental and exploration examples

An elevated concentration alone does not identify contamination, a mineral target
or health risk. Review lithology, sediment fraction, hydrology, land use and human
activity in the appropriate time frame, alongside analytical methods and QA.
These layers are not automatically supplied by the current examples. Unknown
sampling dates, dilution, grain size and QA remain unknown.

Semiquantitative method flags are retained. Uncensored coding is not a claim of
exact measurement precision; reporting-class ranges are not reconstructed from
a method label. Samples sharing coordinates are not automatically certified as
replicates or independent physical sites.

[Return to the project](../README.md) · [Methods](methods.md) · [Sediment examples](examples/README.md)
