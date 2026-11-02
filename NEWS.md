# geoCensor 1.1.0

Released November 1, 2026. This version focuses on USGS National Geochemical
Database workflows.

## Censoring and data integrity

* Parse left and right censoring consistently across estimation and plotting,
  including the USGS `N` qualifier. Preserve observed bounds rather than
  substituting concentrations in censored-model workflows.
* Retain raw values and provenance when recovering rounded upper-limit suffixes
  from the first compatible analysis in a best-value summary.
* Record exclusion reasons for invalid values/limits, conflicting or unknown
  qualifiers, and unsupported concentration units. Normalize equivalent ppm,
  mg/kg, ppb, ug/kg, micro-symbol and percent units.
* Preserve leading zeroes in sample identifiers. Reject duplicate/missing
  metadata keys and prevent joins from multiplying observations. Report
  unmatched and repeated chemistry identifiers.
* Reject duplicate sample/element pairs when widening rather than averaging.
  **Compatibility change:** callers must explicitly select an appropriate
  analysis before creating a wide metals table.

## Estimation and diagnostics

* Support uncensored and right-censored samples in the lognormal estimator.
  Fits require six usable records and two exact detections; count minima are
  computational guards, not evidence of reliable inference.
* Reject fitting warnings and invalid estimates; exclude failed distribution
  fits from AIC rankings.
* Bring `fit_lognorm_with_quants()` under the same interval and numerical fit
  guards as `fit_cen_lognorm()`. Return requested custom quantiles and always
  include the model median. Concentration summaries remain in mg/kg.
  **Compatibility change:** insufficient/failed summaries now have explicit
  `insufficient_data`/`failed_fit` methods instead of an unqualified `MLE` label.
  The successful default output names/order are preserved.
* Add configurable information-review flags without deleting estimates or
  treating their thresholds as validated reliability criteria.
* Add `censored_cdf_bounds()` and `plot_censored_cdf_bounds()` for empirical
  distribution bounds using exact and censored observations, optionally weighted.
  These describe the supplied observations; they are not population confidence
  bands, fitted distributions or significance tests.

## Clustering and sample support

* Add cell-declustering weights averaged over deterministic grid placements and
  a cell-size sensitivity analysis. Report raw exact fractions and the shares
  of weight carried by exact, left-censored and right-censored observations.
* Omit routine confidence limits from declustering-weighted fits. Weight
  effective sample size does not measure independent spatial replication.
* Add `sampling_support_audit()` for metadata coverage, support combinations
  and exact-coordinate review flags. Shared coordinates do not establish
  common physical sites or sample support.
* Retain `spatial_block_bootstrap()` as an explicitly experimental, opt-in
  sensitivity procedure with failure diagnostics and local random-seed
  restoration. Simulations exposed undercoverage and preferential-sampling
  failures; it is not automatically applied to NGDB analyses.
* Keep the covariance-estimation, spatial HAC/calibration and site-block
  spatial-order research outside the installed API. No routine spatial
  confidence procedure or automatic real-NGDB spatial test is introduced.

## Workflows and packaging

* Add `analyze_elements()` and installed command-line examples for single- and
  all-element analyses, with element/domain validation and optional declustering.
* Add a self-contained USGS workflow vignette and cross-platform GitHub Actions
  checks. Full local data runs remain research artifacts, not bundled data.
* Remove active PFAS integration and exclude sibling packages, local research,
  papers and generated outputs from the installed package.
* Correct factor-coordinate parsing in spatial coverage diagnostics and exclude
  blank state labels while preserving numeric coordinate precision.

Ordinary model confidence limits assume independent observations and an
appropriate lognormal model. Spatial confidence limits remain experimental.
