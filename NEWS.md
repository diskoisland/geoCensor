# geoCensor development 1.1.0.9004

* Reject zero-scale likelihood boundaries explicitly in lognormal estimates,
  quantile summaries and fitted plots. When every usable interval includes the
  same uncensored-coded value, numerical roundoff can otherwise produce a
  spurious positive scale and unstable success/failure status. Retain records,
  counts and method review while leaving fitted values unavailable.
* Exclude the same boundary from estimated-scale distribution comparisons with
  an explicit message. Fixed-scale exponential/Rayleigh fits and empirical
  censoring bounds remain available under their existing rules. No precision
  tolerance, reporting interval or concentration is inferred.

# geoCensor development 1.1.0.9003

* Add optional `analyze_elements(..., audit_method_overlap = TRUE)` to append
  element/domain-labelled availability, method coverage and method-pair overlap
  tables. Requires projected coordinates and cell sizes; uses the supplied
  n_offsets. Existing estimates, weights and six/seven result components remain
  unchanged. Can be combined with method stratification independently.
* Add `--method-overlap` to the all-element command, exporting three CSV tables
  and RDS components when requested. The command now requires installed version
  1.1.0.9003 or newer and rejects missing spatial arguments or duplicate switches.
  Missing spatial comparisons remain NA; cell overlap is not paired-sample or
  analytical-agreement evidence, and placement ranges are not confidence limits.

# geoCensor development 1.1.0.9002

* Add `method_spatial_overlap()` to audit occupied/shared grid cells across
  analytical-method labels and spatial scales using a common grid with shifted
  placements. Return record eligibility, method coverage and pairwise Jaccard
  summaries without fitting a concentration model or supplying p-values.
* Retain missing/unspecified methods, excluded chemistry and unlocated records
  explicitly. Comparisons with unavailable spatial data remain NA; observed
  zero overlap requires data for both methods. Grid-placement ranges are
  sensitivity summaries, not confidence limits or evidence of matched samples.
* Preserve distinct representable integer cell indices in the large-grid string
  fallback used by declustering and overlap audits. Default numeric-to-text
  rounding could previously merge neighboring cells at very large grid indices.

# geoCensor development 1.1.0.9001

* Add optional `analyze_elements(..., stratify_methods = TRUE)` summaries by
  normalized analytical-method code within each element/domain. The six existing
  components and their numerical results stay unchanged; an additional
  `method_estimates` table retains counts, fit outcomes and review information.
* Preserve missing/unspecified method buckets and unusable/insufficient/failed
  groups explicitly. These descriptive comparisons do not estimate analytical
  bias or establish comparability across sampling locations or procedures.
* Add `--stratify-methods` to installed single/all-element commands, exporting
  `method_estimates.csv` and the corresponding RDS component when requested.
  Development workflows now require version 1.1.0.9001 or newer.
* Installed commands require a new or empty output directory. Existing results
  are preserved, preventing stale optional tables or plots from appearing in
  a subsequent run with different options or unavailable fits.

# geoCensor development 1.1.0.9000

* Retain analytical-method review in empirical CDF support summaries, CDF plot
  captions and every candidate distribution-comparison row. Invalid or zero-weight
  records are excluded from CDF used counts. Coincident bounds and AIC rankings
  do not establish measurement precision. Numerical results remain unchanged.
* Plot older saved CDF results with method review explicitly unavailable when
  their support metadata predates these fields; do not reconstruct missing counts.
* Carry method-review counts/status through direct lognormal estimates, fitted
  quantile summaries and state summaries, including insufficient/failed fits.
  Review columns are appended; established numeric columns/order are preserved.
  Custom censoring columns and zero-weight exclusions use the actual fit subset.
* Retain method counts and independence assumptions in lognormal plot captions.
  Reject numerical fitting warnings/errors or invalid parameters before plotting,
  consistently with the numerical estimators.

* Add `flag_analytical_methods()` and propagate method review through long
  observations, batch estimates, declustering summaries and installed CLI exports.
  Documented USGS semiquantitative codes ES_SQ and ES_H2O_P receive review flags;
  missing, unspecified and unreviewed method information remains unresolved/NA.
* Add a `methods` table to `analyze_elements()` with per-code input and unweighted
  fit-subset counts, review status and provenance. Method concerns are separate
  from chemistry validity and count-based review flags. Numerical estimates,
  censoring fields and record inclusion are unchanged; no class bounds are guessed.
* Declustering method counts use the same valid chemistry/coordinate subset as
  each fit, rather than borrowing counts from the unweighted estimate.
* Single-element CLI output now uses the shared public batch workflow and exports
  `methods.csv`; estimate/model tables gain additive element/count/review columns.
* Clarify that n_exact is an uncensored coding count, not verified measurement
  precision. Correct the vignette's soil horizon-zone example to `sample_zone`.
* Keep installed-source fallback dependencies complete for the Utah dust example.

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
