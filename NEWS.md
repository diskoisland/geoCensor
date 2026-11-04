# geoCensor 1.2.0

## Installed workflows and worked examples

* Add an Upper Illinois River Basin vignette and installed case-study script.
  Select all reported species by the project metadata, retain original source
  fields and record ordinals, audit support, compare 5/10/20 km declustering,
  and save results, figures, provenance and unavailable fits. The actual source
  data are downloaded separately; no private research bundle is included.
* Accept the final 1.2.0 version in installed commands, while checking the
  required exported functions and arguments to reject incomplete earlier local
  candidates with the same version label. Restart R after upgrading.
* Spatial intervals remain experimental, and routine declustering-weighted
  confidence limits remain absent. The case study does not validate geological
  background, ore probabilities or contamination attribution.

## Disk exports

* Add optional `export_usgs_chemistry()` and an installed `export_chemistry.R`
  command for batch exports without retaining selected data frames. All original
  text columns and separate source-record ordinals are written to a new directory.
  Source checksums, import provenance and output checksums accompany the export.
* Reserve destinations for cooperating calls, reject existing outputs, and publish
  only after complete parsing and matching before/after source fingerprints.
  Existing readers and analysis defaults are unchanged.

## Species-selective imports

* Add optional `read_usgs_chemistry_select()` to scan all input batches while
  retaining requested species, original text fields and separate source-record
  ordinals. Malformed unselected or late records still abort the read. Parser
  batches are bounded; retained output still grows with selected-record count.
* Add explicit `--chunk-size=10000` to the single-element analysis command,
  with a separate input-record index. Its default whole-file reader is unchanged.

## Lossless CSV import

* Read chemistry CSV columns as text with `read_usgs_chemistry()`, shared by
  `build_ag_from_csv()` and both installed analysis commands. Preserve rare
  qualifiers, upper-limit suffixes, lexical numeric precision and identifier
  formatting. Reject malformed CSV structure and duplicate/empty headers.
* Fix a sparse-column import failure that could drop a rare positive-value G/L
  qualifier and turn a censoring limit into an apparent detection. Source-review
  exports now retain character source fields; in-memory APIs are unchanged.
* Require installed 1.2.0 for the analysis commands. Read command metadata
  as text and explicitly parse only requested projected coordinates; invalid
  coordinate text stops before saving output. Explicit nonfinite coordinates
  remain unlocated. No CRS, units or missing data inferred.

## Saved sampling-support audits

* Add `audit_sampling_file()` and installed `audit_sampling_support.R` for
  reproducible metadata audits and optional position memberships. Export raw
  metadata, review tables, coordinate-conversion provenance, RDS and run/session
  records. Preserve raw text and leading-zero IDs; declare analysis missing markers.
* Parse only explicit numeric coordinate text, retaining nonfinite/missing values
  as unlocated and rejecting invalid text before output creation. Preserve
  nonempty destinations and check input fingerprints before saving. No coordinate
  system, units, dates, physical identity, weights or dilution factors inferred.

## Recorded-position memberships

* Add optional `sampling_support_audit(..., include_memberships = TRUE)` with
  original-row memberships and exact-coordinate groups across analysis domains.
  Preserve raw IDs/fields and unlocated rows. Missing coordinates never form
  one artificial group, and unavailable group sizes remain NA.
* Keep the four default components unchanged. New cohort-specific position IDs
  do not establish physical sites, replicates or independence; no observations
  are collapsed or weighted. Exact numeric grouping avoids text-rounding merges.

## Guarded upper-limit recovery

* Require a known matching analytical method, well-formed first entry and
  compatible usable coding before recovering rounded summary upper limits.
  Missing/UNKNOWN methods, mismatches and conflicting or malformed evidence
  now block recovery. Original qvalue and qualifier are retained and parsed
  according to their own coding; another method's limit is never borrowed.
* Add `summary_recovery_status` to builder observations to distinguish recovery,
  original upper coding, absent evidence and specific blocking reasons. This
  intentional compatibility change can alter previous unguarded recoveries;
  correctly associated recovery and source-precision auditing remain available.
  No dilution correction or analytical-QA certification is inferred.
* Require installed 1.2.0 for both command workflows, including runs without
  source auditing, so new scripts cannot silently use an older unguarded parser.

## Batch source-summary review

* Integrate the source-summary precision audit into the active package. Add
  opt-in `analyze_elements(..., audit_source_precision = TRUE)` row-level review
  and element/domain/status counts. Retain original requested chemistry rows,
  native units and raw IDs; preserve all existing estimation outputs.
* Expose `--audit-source-precision` in both installed analysis commands, exporting
  the two additional tables and RDS components only when requested.
* Treat literal UNKNOWN method labels as unresolved association in the audit,
  including when both reported and summary labels say UNKNOWN. No source-value
  candidate is exposed from unspecified method identity.

## Source-summary precision audit

* Add `audit_summary_precision()` to expose first-entry source numeric differences and upper-limit evidence only when method association and censoring checks permit. Preserve every input row and value; candidates are review evidence and never replace estimator inputs.

## Analytical methods and geographic support

* Add `flag_analytical_methods()` and retain review counts/status through parsed
  observations, estimates, quantiles, plots, empirical CDF support and distribution
  diagnostics. Documented USGS ES_SQ and ES_H2O_P codes are flagged; missing,
  unspecified and other supplied codes remain unresolved or unreviewed. Raw
  labels and chemistry are retained; no reporting-class bounds are inferred.
* Add a `methods` component to `analyze_elements()` and optional
  `stratify_methods = TRUE` summaries within each element/domain. Missing and
  UNKNOWN labels remain separate unresolved buckets, including insufficient
  and failed groups. Method contrasts are descriptive, not calibrations.
* Add `method_spatial_overlap()` and optional
  `analyze_elements(..., audit_method_overlap = TRUE)` with projected geometry.
  Return availability, method coverage and shared/union-cell/Jaccard summaries
  across fixed grid placements. Missing comparisons remain NA; observed zero
  overlap requires spatial data for both labels. Shared cells do not establish
  paired samples, analytical agreement or comparable geological support.
* Expose method stratification in the single/all-element commands and geographic
  auditing in the all-element command. Both commands require geoCensor 1.2.0 or
  newer. Optional output tables appear only when requested; existing output
  directories are preserved rather than reused or overwritten.

## Numerical safety and interpretation

* Reject zero-scale likelihood boundaries explicitly. When all usable intervals
  include a single uncensored-coded value, the model has no finite positive-scale
  maximum. Estimates/quantiles remain missing and fitted curves are unavailable,
  even if numerical roundoff could otherwise report a tiny positive scale.
  Counts, raw records and method review are retained.
* Exclude these boundaries from estimated-scale distribution comparisons with
  an explicit message. Fixed-scale exponential/Rayleigh models and empirical
  censoring bounds retain their ordinary rules. No precision tolerance is used.
* Preserve distinct representable integer grid cells at very large indices;
  default numeric-to-text rounding could previously combine neighboring cells.
* Preserve method-review cautions in saved plots. Older saved empirical-CDF
  objects report method review unavailable rather than reconstructing counts.
* Extend the runnable USGS vignette with available, all-left, constant-coded and
  missing-method examples. Method-specific declustering remains separate research
  outside the installed API. Weighted routine confidence limits stay missing;
  spatial confidence limits remain experimental.

## Compatibility

Existing numerical results are preserved apart from the explicit fitting-boundary
and large-grid corrections. Review columns and the default `methods` component
are additive; code should select named columns/components. Optional method and
geographic tables do not select records, methods or a preferred spatial scale.
“Uncensored” and “exact” refer to input coding, not verified measurement precision.
An available fit or cleared count-review flag does not establish model adequacy,
independence, a reporting design or a geological population.

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
