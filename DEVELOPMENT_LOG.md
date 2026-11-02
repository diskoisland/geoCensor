# geoCensor development progress

## Status — November 1, 2026

The package source version is **1.1.0**, consistently recorded in
`DESCRIPTION`, the README citation and `CITATION.cff`. Release record and checked
assets: [v1.1.0](https://github.com/diskoisland/geoCensor/releases/tag/v1.1.0).
Historical entries below retain the earlier 0.1.0 development version and the
prepublication candidate status. The starting commit was `911f8aa` (Add censored
lognormal fit plots). Publication is limited to an explicit file allowlist.

The candidate passed 77 tests / 410 expectations without failures, errors,
warnings or skips, plus installed single/all-element examples. All five remote
jobs passed on both candidate and documentation heads. Final release wording
received a fresh clean source build and R CMD check (Status: OK), including
vignette creation/rebuilding and package tests. The R/help/test/example/vignette
sources remain unchanged from the validated candidate. Final exact-commit
publication checks are specified in [the release validation record](RELEASE_CHECKLIST.md)
and retained with the GitHub release assets.

Completed package work includes sampling-support/method audits, multi-element
USGS analysis, information review flags and a worked vignette. Spatial research
has progressed through Gaussian calibration, weighted censored-score transfer,
and independently changing survey layouts. The latter exposed heavy-censoring
bias and undercoverage. The current focused work is a spatial censored-likelihood
reference, descriptive censoring-weight support, and empirical CDF censoring
bounds, with numerical and package
validation recorded in the dated entries below. The candidate source is now 1.1.0; historical research entries preserve their
original version/status records.

Latest package validation: the censoring-weight shares and new empirical CDF
bounds pass the full local test suite and R CMD check (Status: OK), including
vignette rebuilding. The
all-USGS support audit reconciles 22,431 rows across all 80 soil/81 sediment
species. Experimental likelihood work has passed numerical gates but still
shows sensitivity when scale is estimated under heavy censoring; routine
spatial confidence limits remain withheld.

A fresh estimated-scale pilot now covers 16 new survey layouts and 64 paired-
censoring analyses. Numerical checks were stable, while sparse detections left
some estimates sensitive. The worst two cases fail the existing public exact-
count rule, and the secondary count-eligible results are reported explicitly.
This is point-estimation research with known covariance shape, not calibrated
spatial inference; see the final dated entry below.

Selected location/scale profiles now show broad, asymmetric likelihood support
in two count-eligible heavy-censoring cases, compared with sharper moderate-
censoring partners. All 324 profiled maxima are finite and the independent
numerical audits are stable. The existing few-exact review warnings still apply;
no new confidence procedure or exclusion rule has been introduced.

The next fresh pilot now estimates location, marginal SD, exponential range and
nugget fraction from observations. Its original sixteen attempts give thirteen
interior and three search-boundary outcomes; the sole flagged-interior case
moves to a boundary on a separate higher-precision check. Twelve primary cases
remain unflagged interior outcomes. Unresolved accepted values remain missing,
and a separate noise-only endpoint check records unidentified range as NA.
Confidence limits remain experimental; no public covariance estimator was added.

The explicit-endpoint revision has now completed another fresh study: four
layouts, twelve fields and twenty-four analyses. Sixteen comparisons separate
numerically and eight remain unresolved; overall range is missing in fourteen
cases (six independent-noise results plus eight unresolved). All endpoint branch
ranges are NA. Positive range also appears on seven endpoint-generated cases,
so numerical likelihood preference cannot be treated as detected spatial signal.
That statistical evidence check has now passed a separate restricted-null fresh
validation, described below. Real-USGS eligibility remains the next gate;
confidence limits stay experimental.

The real-data eligibility audit now covers both export dictionaries/schemas
and all readable soil metadata. It identifies substantial repeated-coordinate
support ambiguity and incomplete reporting-design documentation. Soil counts
were independently reconciled; chemistry and sediment record reads timed out,
so their newly audited contents remain unavailable. No real-data p-value or
count of eligible domains is asserted; see the latest dated audit entry.

This log records implementation progress and validation. [NEWS.md](NEWS.md)
contains the user-facing change list. Append dated entries as work continues.

## Completed in this work session

### Package scope and builds

- Focused geoCensor on USGS geochemical data by removing the PFAS integration
  from package functions, exports, help files, examples, and tests.
- Removed PFAS-only `readxl` and `NADA2` dependencies. The separate `pfascensor`
  folder and data were preserved outside the package build.
- Extended `.Rbuildignore` to exclude sibling packages, local analyses,
  generated documents, outputs, check folders, and temporary files.

### Input integrity

- Added explicit `validation_status` reasons for invalid or missing values,
  nonpositive detections, invalid limits, conflicting censoring encodings,
  unknown qualifiers, and unsupported units.
- Standardized concentration conversions for ppm/mg/kg, ppb/ug/kg (including
  micro-symbol variants), and percent units.
- Preserved leading zeroes in chemistry CSV sample IDs.
- Rejected duplicate or blank metadata IDs before joins; warned about unmatched
  and repeated chemistry IDs, retained those chemistry rows, and added a
  `join_audit` attribute with counts.
- Replaced silent averaging of duplicate sample/element results in wide tables
  with an error requiring explicit result selection. Added validation summaries.

### Censoring and model reliability

- Enabled lognormal estimates for uncensored and right-censored-only samples,
  subject to the existing minimum of six usable rows and two exact detections.
- Unified interval construction across estimation, diagnostics, and plotting.
- Prevented invalid censoring limits from turning flagged rows into detections.
- Returned missing estimates for fitting errors, warnings, or invalid parameters;
  excluded unreliable fits from AIC ranking and Akaike weights.
- Recognized the USGS `N` qualifier as left censoring.
- Full-data inspection revealed rounded CSV `qvalue` fields that hid upper-limit
  suffixes. Added recovery of explicit `.1111`/`.11111` markers from the first
  best-ranked `summary` analysis only when it agrees with rounded `qvalue`.
  Retained `qvalue_raw` and recorded the source in `censoring_source`.

### Spatial clustering

- Added `cell_declustering_weights()` for projected coordinates. Samples share
  their occupied cell's total influence; weights are averaged across 100
  deterministic grid placements by default.
- Added optional observation weights to `fit_cen_lognorm()`, retaining censoring
  intervals and omitting ordinary confidence limits for weighted results.
- Added `declustering_sensitivity()` to compare several cell sizes against an
  unweighted baseline on the same usable chemistry/coordinate subset.
- Added a weight effective-sample-size diagnostic, explicitly distinguished
  from a count of spatially independent observations.
- Reduced large-run allocation costs with numeric cell identifiers when exact
  integer representation is possible, retaining string identifiers otherwise.

### Multi-element workflows

- Added `analyze_elements()` to analyze all reported species or a selected list,
  separately within optional state, medium, or geological-domain groups.
- Included estimates, model comparisons, validation counts, optional spatial
  sensitivity, and join audits. Insufficient-data species remain in results.
- Recalculated declustering weights for each element/domain's usable data.
- Added installed examples for single-element and all-element command-line runs.
- Added checkpointed full-data and result-summary scripts under
  `scripts/analysis/`, outside the installable package.

## Validation completed

- The full local test suite passed, including new regression tests for input
  integrity, model failures, declustering, and multi-element/domain isolation.
- End-to-end single-element and batch workflows passed on synthetic data,
  including unsupported units, duplicate samples, missing species, and grouped
  spatial runs.
- A clean source build and `R CMD check --no-manual --no-build-vignettes`
  completed with **Status: OK** after the full-data import correction.
- The check used installed dependencies; remote CRAN/Bioconductor indexes were
  unavailable. PDF manual generation and a separate vignette build were not
  validated by that command.
- Full-data output checks reconciled input, validation, and grouped row counts;
  verified three spatial scales per state/species; and verified missing weighted
  confidence limits as intended.

## Full USGS run completed

Inputs were the local `ngdbsoil` and `ngdbsed` best-value and metadata CSVs.
Source file contents were not modified. Only samples explicitly marked United
States were included; reported state/territory groups were retained. Soil and
sediment were processed separately.

| Measure | Soil | Sediment |
| --- | ---: | ---: |
| Original chemistry records | 3,680,506 | 10,567,162 |
| U.S. chemistry records analyzed | 3,616,489 | 10,175,588 |
| U.S. metadata sample records | 102,035 | 322,988 |
| Reported species analyzed | 80 | 81 |
| Finite pooled estimates | 76 | 79 |
| Finite state/territory estimates | 2,836 | 3,113 |
| Finite weighted spatial estimates, across three scales | 8,511 | 9,327 |
| Recovered upper-censoring markers | 5,040 | 22,296 |
| Chemistry records excluded by validation | 876 | 10,202 |

Reported species include compounds and ancillary sample-weight fields, not just
element symbols. Sample weights reported in grams are not concentration data.
Soil osmium and rhenium also lacked enough exact detections for pooled fits.

Spatial sensitivity used 5, 10, and 20 km cells with 100 grid placements per size.
Each state/territory used a local WGS84 Lambert azimuthal equal-area projection;
definitions were saved alongside the outputs. The scales are exploratory, not
geologically optimized sampling spacings.

Results are under `output/usgs_all_elements_2026-10-31/`:

- [Run summary](output/usgs_all_elements_2026-10-31/README.md)
- [Run counts](output/usgs_all_elements_2026-10-31/run_overview.csv)
- [Combined spatial comparison](output/usgs_all_elements_2026-10-31/declustering_comparison.csv)
- `soil/` and `sed/`: pooled and grouped estimates, distribution comparisons,
  validation tables, spatial results, join audits, projections, and checkpoints.

These generated results are local and excluded from package builds. An
interrupted preliminary soil run was archived separately under
`output/usgs_all_elements_2026-10-31_preliminary_before_censoring_QC/`;
use the completed, corrected run above.

Reproduction scripts:

- [Full-data run](scripts/analysis/run_all_usgs_elements.R)
- [Output checks and summary](scripts/analysis/summarize_all_usgs_elements.R)

The runner resumes its checkpoints. Use a new output directory for a changed
model, input dataset, or analysis configuration so stale checkpoints are not
reused.

## Spatial uncertainty pilot — October 31, 2026

Implemented `spatial_block_bootstrap()` as an explicitly experimental API.
It resamples occupied nonoverlapping blocks, recalculates declustering weights,
preserves block-draw multiplicity, refits the censored model, and returns
percentile limits with replicate and block diagnostics. It suppresses limits
for too few blocks, unavailable point fits, or excessive failed replicates.
Seeded calls restore the caller's random-number state. Existing weighted
estimates and full-data outputs retain their previous treatment of uncertainty.

Software regression tests and the full local suite passed. A fresh source build
and `R CMD check --no-manual --no-build-vignettes` finished with **Status: OK**.

Completed 100 simulated datasets in each of five prespecified settings, with
199 bootstrap replicates per dataset (99,500 successful bootstrap fits).
Both left and right censoring were included. Nominal interval coverage was 95%.

| Simulation setting | Observed coverage |
| --- | ---: |
| Independent sites, regular density | 88% |
| Independent sites, uneven density | 91% |
| Independent sites, grade-dependent density | 94% |
| Correlated sites, smaller blocks | 73% |
| Correlated sites, larger blocks | 82% |

**Decision: the prototype is not validated for routine NGDB confidence intervals.**
Successful model convergence did not establish calibrated interval coverage.
Increasing block width helped the correlated setting but did not reach nominal
coverage. The study used 64 observed sites on a grid and did not simulate missing
regions or preferential omission of entire sites. It is a pilot, not a complete
operating-characteristic study for irregular NGDB sampling.

Results and reproducible scripts:

- [Pilot report](output/spatial_bootstrap_validation_2026-10-31/README.md)
- [Coverage summary and Monte Carlo intervals](output/spatial_bootstrap_validation_2026-10-31/coverage_summary.csv)
- [Coverage figure](output/spatial_bootstrap_validation_2026-10-31/coverage.png)
- [Simulation runner](scripts/analysis/validate_spatial_bootstrap.R)
- [Report and output checks](scripts/analysis/report_spatial_bootstrap.R)

Next statistical work is to compare a design-appropriate spatial bootstrap or
spatial model, test interval calibration, and expand the simulation designs.
No bootstrap intervals were added to the full USGS analysis.

## Bootstrap benchmark and diagnosis — October 31, 2026

Completed the next diagnostic step before changing production uncertainty methods.
A new regression test confirms exact agreement between uncensored package
bootstrap replicates and resampled site log means. The full test suite and a
fresh local package check passed (Status: OK).

The prespecified study used 2,000 development plus 3,000 independently seeded
confirmation datasets for the exact independent Gaussian site-mean benchmark
(999 resamples each). Censored diagnostics used 200 development plus 200
confirmation datasets per setting, with 199 actual package refits each.
Candidate formulas and scenarios were fixed before confirmation; no method was
selected or retuned on the development results.

Primary target: the superpopulation geometric mean exp(E(log concentration)) =
50. The realized equal-site domain target was evaluated separately and never
substituted to improve apparent coverage. Neither target establishes a mean for
unobserved terrain.

| Independent confirmation result | Coverage of nominal 95% interval |
| --- | ---: |
| Uncensored analytic t benchmark (3,000 datasets) | 95.13% |
| Uncensored studentized bootstrap (3,000 datasets) | 95.00% |
| Uncensored current percentile method (3,000 datasets) | 94.57% |
| Censored independent sites, regular density (200 datasets) | 96.0% |
| Censored independent sites, uneven density (200 datasets) | 96.5% |
| Correlated sites, smaller blocks (200 datasets) | 77.0% |
| Correlated sites, larger blocks (200 datasets) | 85.0% |

The original 88% independent-site pilot result did not persist in the larger
simple benchmark. Censoring and uneven density alone did not reproduce the
major shortfall. Under correlation, reflected basic and centered normal log-scale
intervals also undercovered. Point-estimate log bias was approximately 0.001;
bootstrap standard errors were only 0.63 and 0.76 times the empirical SD.

**Diagnosis within these scenarios:** dependence across blocks is the principal
unresolved source of underestimated uncertainty. Changing the interval formula
alone is insufficient. This does not prove general reliability for irregular or
preferentially incomplete sampling. All 64 sites were observed; blocks and
correlation scales were fixed, and censored confirmation coverage has Monte
Carlo uncertainty reported in the output tables.

Production intervals remain experimental. Next compare a design-appropriate
resampling approach or spatial covariance model, validate uncertainty under
irregular locations and preferential omission of sites, and define the intended
geological population before applying routine intervals to NGDB results.

- [Diagnostic report](output/spatial_bootstrap_diagnostics_2026-10-31/README.md)
- [Detailed coverage, bias, and SE checks](output/spatial_bootstrap_diagnostics_2026-10-31/summary.csv)
- [Confirmation figure](output/spatial_bootstrap_diagnostics_2026-10-31/confirmation_coverage.png)
- [Study runner](scripts/analysis/diagnose_spatial_bootstrap.R)
- [Report and output checks](scripts/analysis/report_spatial_bootstrap_diagnostics.R)

## Remaining work toward V1.1

- Validate spatial uncertainty beyond the initial block-bootstrap pilot before
  adding routine weighted confidence intervals. The current prototype undercovers.
- Evaluate cell sizes and geological domains, and distinguish true coverage of
  a defined area from spread among sampled locations. Current declustering does
  not estimate concentrations in unsampled terrain.
- Review method, horizon, depth, and size-fraction differences before interpreting
  pooled populations. Distinct sample IDs at the same location remain separate.
- Review low-information estimates using detection counts and cell-size
  sensitivity. Minimum row counts alone do not establish reliable inference.
- Review and push the prepared cross-platform checks, then verify the new USGS
  vignette and tests on those remote platforms. Local builds already pass.
- Review and commit the implementation, then synchronize version and citation
  metadata, finalize release notes, and create the intended V1.1 release/tag.

Kriging, entropy-based sampling allocation, and automatic exploration-budget
recommendations have not been implemented. Ordinary pooled confidence intervals
remain model-based and do not account for spatial dependence.

## Authorized autonomous continuation — October 31, 2026

The user authorized continued package development while resting, with progress
recorded here. Work will preserve source datasets and existing changes; no
release will be published in this phase.

Planned sequence: audit real sampling-support and analytical-method differences;
add a reusable, descriptive comparability audit; extend uncertainty validation
to irregular sites and genuinely preferential omission of sites; improve package
checks and the worked example where justified. Uncertainty methods remain
experimental unless the evidence supports a narrower validated use.

### Autonomous phase: support audit and reproducible workflow

Implemented and tested `sampling_support_audit()` for field missingness,
recorded support combinations, and exact-coordinate reuse. It reports differing
recorded support separately from missing information and does not merge or
remove samples. Explicit field/ID validation prevents accidental inflation.
The use of dplyr::pick is reflected in dplyr >= 1.1.0 in DESCRIPTION.

Completed the descriptive audit of all 102,035 U.S. soil metadata records and
322,988 sediment metadata records. Counts reconcile with combinations and
coordinate summaries. Soil has 12,559 reused coordinate pairs containing 55,226
sample records, with differing recorded support at 5,504 locations. Sediment
has 26,242 reused pairs containing 91,808 records, with differing support at
5,251 locations. Exact coordinate reuse is not evidence of duplicate samples.

Depth is recorded for 12.64% of soil samples and 4.02% of sediment samples;
soil horizon for 17.04%. Mesh size is recorded for 24.93% and 38.43%, respectively.
Multiple analytical method codes occur in 2,421 of 3,697 soil species/state groups
and 2,864 of 3,780 sediment groups. No automatic support/method exclusion was
introduced; missing metadata does not establish equivalence or justify a
convenient imputed sampling support.

- [Support and method audit report](output/usgs_sampling_support_audit_2026-10-31/README.md)
- [Audit totals](output/usgs_sampling_support_audit_2026-10-31/overview.csv)
- [Audit runner](scripts/analysis/audit_usgs_sampling_support.R)
- [Report and reconciliation](scripts/analysis/report_sampling_support_audit.R)

Prepared a local GitHub Actions workflow for release R on Linux/macOS/Windows
plus old-release and development R on Linux. It uses the official r-lib/actions
v2 pattern, with read-only repository permissions. It has not been pushed or
run on GitHub; local YAML parsing is not cross-platform execution validation.

Added a self-contained synthetic USGS workflow vignette covering support review,
censoring recovery, joins, domain selection, all-element estimates, spatial
sensitivity, experimental interval limitations, and reproducible exports. Its
chunks ran and the HTML rendered successfully against current package source.
The first direct render loaded the older installed package; source-loaded
rendering corrected that environment mismatch. Pandoc emitted an upstream
highlight-style deprecation message without preventing output.

An independently seeded irregular-location/preferential-site-omission study is
running in a fresh output directory. No production spatial intervals have been
promoted or assigned to the full USGS outputs.

### Autonomous phase: irregular selection, information flags, and build checks

Completed the irregular-site study: 100 development plus 100 independently seeded
confirmation datasets per each of five prespecified settings, with 199 package
resamples each (199,000 successful refits). Only 64 of 144 jittered candidate
sites were observed. Counts within observed sites varied independently of grade.

Confirmation coverage of nominal 95% intervals was 96% for independent uniformly
selected sites and 87% for correlated uniformly selected sites. Under the
specified high-concentration site-selection mechanism, coverage was 0% (100
confirmation datasets; Monte Carlo upper 95% bound 3.62%) and the mean estimate
was 75.74 for a population target of 50. The mean complete-data observed-site
log bias relative to the full realized domain was about 0.401, showing that the
selection error remains even without censoring. Shifting block origins did not
repair it. These are controlled synthetic results, not estimates of true NGDB
selection bias.

The uniform correlated setting had 99% coverage for the separately defined
realized-domain target but 87% for the population target. This reinforces the
need to define the estimand rather than changing targets to improve an interval's
reported performance. No new routine uncertainty method is justified by these
results; experimental labels and missing routine spatial intervals are retained.

- [Irregular selection report](output/irregular_spatial_validation_2026-10-31/README.md)
- [Study settings, bias, and coverage](output/irregular_spatial_validation_2026-10-31/summary.csv)
- [Simulation runner](scripts/analysis/validate_irregular_spatial_sampling.R)
- [Report and checks](scripts/analysis/report_irregular_spatial_sampling.R)

Added `flag_estimate_information()` and optional batch review thresholds. Defaults
review fewer than 10 exact observations or less than 5% exact observations;
these are configurable count conventions, not validated inferential cutoffs.
Right-censored values are excluded from the exact count, not relabeled as
nondetections. No estimates or records are removed or recalculated.

Annotated the completed full-data outputs in a separate folder, checking every
geometric mean against its original value. Count flags request review for 5
available pooled soil estimates and 9 sediment estimates, and for 452 soil and
384 sediment state/species estimates. Spatial rows repeat elements/domains across
cell sizes and must not be counted as independent analyses.

- [Information review notes](output/usgs_information_review_2026-10-31/README.md)
- [Review counts](output/usgs_information_review_2026-10-31/review_counts.csv)
- [Annotation runner](scripts/analysis/annotate_full_run_information.R)

Generated outputs, check directories, and package archives are now ignored by Git
by default; files were retained on disk. Source/tests/documentation remain
reviewable. Unrelated working documents and sibling projects were not removed.

A clean source build and `R CMD check --no-manual` passed with the new workflow
vignette built and rebuilt. Remote repository indexes were unavailable during
local dependency checking; installed dependencies were used. The GitHub workflow
has not run remotely. PDF manual generation has not been checked. A final check
will follow the information-flag integration and documentation updates.

### Autonomous phase: final verified checkpoint

Final clean source build and `R CMD check --no-manual` completed with **Status:
OK**, including the full test suite and vignette building/rebuilding after the
information-flag integration. Also fixed spatial coverage parsing of factor
coordinate labels so factors cannot be mistaken for numeric level codes; blank
state labels are excluded. The regression test passed and numeric coordinate
precision is preserved for already numeric inputs.

- [Tested development snapshot and review notes](output/autonomous_build_validation_2026-10-31/README.md)
- [Saved check log](output/autonomous_build_validation_2026-10-31/00check.log)
- [Built workflow example](output/autonomous_build_validation_2026-10-31/usgs-workflow.html)

The snapshot archive and its SHA-256 identify the checked source. The local
check used installed dependencies because remote indexes were unavailable.
PDF manual generation and actual remote CI execution remain unverified. No
source files were staged, committed, pushed, or published during this pass.
Generated outputs and source datasets remain local and preserved.

The next scientifically justified phase is to define a comparable reference
population using the support/method audit and then compare a design-appropriate
spatial/selection model or resampling method under irregular and informative
sampling. A covariance-only correction cannot by itself repair omitted-site
selection bias. Avoid introducing routine geographical-domain confidence
intervals simply because a toy simulation can be tuned to nominal coverage.

Current deliverables are reviewable without deleting or replacing earlier
results. The existing full-data estimates remain exploratory, and all new
information annotations preserve their original geometric means.

## Conley spatial HAC investigation — October 31, 2026

The user explicitly directed attention to Conley/HAC and supplied
[SpatialInference's bandwidth vignette](https://cran.r-universe.dev/articles/SpatialInference/spatial_HAC.html).
Read the rendered article, official reference manual, Conley (1999), and the
sandwich/survival score documentation. SpatialInference's high-level wrappers
are documented for linear/felm models, not a direct replacement for the weighted
interval-censored survreg estimator. The package is not installed locally and
its C++ output has not yet been used as a numerical comparison.

The useful transferable ideas are distance-based covariance of estimating
functions, kernel/cutoff sensitivity, and covariogram-based range diagnostics.
A very wide cutoff is not automatically conservative. In a score sandwich,
a near-all-ones kernel can collapse the meat toward the outer product of the
summed fitted scores, which are near zero. This behavior was verified directly.
Bandwidth diagnostics for the censored estimator should be based on the relevant
weighted score/influence contributions, not OLS residuals of substituted values.

Built a research-only reference under scripts/analysis/conley_survreg_reference.R.
It uses the full joint location/log-scale weighted likelihood scores and inverse
weighted observed information. Co-located score vectors are summed before S' K S.
The dense reference has a 2,000-location guard and is not a national-data engine.
No public package interval API or full-data result was changed.

Reference checks passed: numerical likelihood derivatives for both parameters,
bread agreement with survreg, isolated-location HC0 reduction, arbitrary weight
scale invariance, the uncensored location-mean variance identity, and wide-cutoff
SE shrinkage. Scripts and saved results:

- [Score/reference verification](scripts/analysis/verify_conley_score_reference.R)
- [Reference calculation](scripts/analysis/conley_survreg_reference.R)
- [Checks](output/conley_reference_checks_2026-10-31/checks.csv)
- [Local score-extractor comparison](output/conley_reference_checks_2026-10-31/score_extractor_comparison.csv)

A concrete compatibility issue was reproduced locally with survival 3.8.6 and
sandwich 3.1.3: sandwich::estfun() applies case weights to the location score but
its scale column did not match the numerical derivative of the weighted
likelihood in this test. The direct survival matrix derivatives, multiplied by
weights for both parameters, agreed to about 1e-10. This is a version- and
case-specific verified observation, not a claim that every sandwich method is
incorrect. No upstream bug report or external message was sent.

A prespecified study is now running with independent development/confirmation
seeds: 200 datasets per phase/setting; 64 or 256 observed sites; independent,
correlated, and preferential-omission settings; fixed bandwidths and product
Bartlett/radial Bartlett/Epanechnikov kernels. Negative or non-PSD covariances
are recorded rather than silently repaired. No kernel or cutoff will be selected
on coverage results from this study alone.

### Conley phase: completed reference study and external numerical check

Completed the prespecified score-HAC study: 200 development plus 200 confirmation
datasets in each of four settings (1,600 fits), six fixed bandwidths and three
kernels (28,800 covariance evaluations). Numerical/reference checks passed.
No kernel or bandwidth was chosen as a production rule from these results.

In confirmation, product Bartlett coverage in correlated 64-site data was
79.5% at bandwidth 4 and 57.5% at bandwidth 16. With 256 correlated sites,
coverage was 89.5% at bandwidths 4 and 8. Epanechnikov at bandwidth 4 gave
92.5% in the 256-site setting; some wider Epanechnikov covariance matrices were
not positive semidefinite and were withheld, not clipped. Preferential omission
retained large point bias and zero coverage across these candidates.

Installed SpatialInference 0.1.0 (Remote SHA
5bc0beeb63d19fc369cae304021b2a6b1440aca9) into an isolated temporary library for
reference comparisons, without changing geoCensor dependencies or the user's
normal library. Its dense XeeXhC core agreed with S' K S to about 1e-15 across
Bartlett, Epanechnikov, uniform, Gaussian, Parzen, and biweight kernels in the
specified four-location/two-column fixture. This validates computational algebra,
not coverage or direct compatibility of its linear wrappers with survreg.

The documented XeeXhC_Lg call produced a matrix-dimension error for that fixture
in this exact installed version. A minimal reproducible comparison is saved;
no global claim about other versions was made and no upstream message was sent.
The dense engine is therefore a useful numerical reference, while a scalable
adapter needs separate verification before any national-data use.

- [Conley/HAC study report](output/conley_hac_validation_2026-10-31/README.md)
- [Coverage and covariance failures](output/conley_hac_validation_2026-10-31/summary.csv)
- [External core comparison](output/conley_reference_checks_2026-10-31/spatialinference_core_comparison.csv)
- [External comparison script](scripts/analysis/compare_spatialinference_core.R)
- [Study runner](scripts/analysis/validate_conley_reference.R)
- [Report and output checks](scripts/analysis/report_conley_reference.R)

Next justified work: investigate influence-based bandwidth diagnostics and
finite-sample/fixed-bandwidth calibration with a separately prespecified study;
compare the package's reusable computational core after a scalable adapter is
verified. Keep the covariance problem distinct from informative site selection.
The same M-estimator consistency requirement discussed in the sandwich reference
still applies: a covariance correction does not make a biased estimating equation
unbiased. Routine full-data spatial intervals remain unchanged and experimental.

## Literature handoff and Gaussian-reference calibration — October 31, 2026

Searched primary publisher, author, and preprint sources and prepared
[PAPER_FETCH_LIST.md](PAPER_FETCH_LIST.md): 13 prioritized papers, DOI/open-PDF
links, access observations, and explicit implementation questions. The first
four published fetch priorities are Bester et al. (2016) on spatial fixed-b
inference; Conley (1999); Ordoñez et al. (2018) on censored spatial likelihood;
and Qiang & Craigmile (2026) on preferential sampling. Open alternatives are
linked for the other priorities. Abstract/metadata screening is distinguished
from full-text reading; no claim of full review is made. This handoff is excluded
from package builds.

Moved to a new research experiment rather than waiting for papers. It isolates
critical-value calibration from HAC covariance estimation. Gaussian site-mean
reference distributions use fixed assumed exponential covariance classes, 10,000
calibration draws and 5,000 independently seeded validation draws per design/range.
These assumptions are explicit and not claimed known for NGDB.

The transfer experiment fitted 4,800 censored datasets (300 development plus 500
independent confirmation datasets for each of six site-count/range settings),
evaluating 76,800 intervals. No method or assumed range was selected or retuned
after confirmation. Product Bartlett bandwidths were fixed at 2, 4, 8, and 16.
This is an illustrative Gaussian-reference experiment, not an implementation
of published Müller–Watson or fixed-b methods.

At the prespecified bandwidth 4 and correctly matched reference decay parameter
0.8, censored confirmation coverage rose from 83.8%/86.0% under 1.96 to
96.2%/96.4% under reference-calibrated critical values for 64/256 sites. When
actual dependence had decay parameter 2 but the assumed reference remained 0.8,
coverage fell to 83.8%/84.8%. Taking the maximum critical over the finite assumed
range grid produced 96.4%/95.2% there, at the cost of conservative or 100% coverage
under weaker dependence. It does not establish protection outside the assumed
Gaussian class, under irregular/informative selection, or for unknown support.
Routine package intervals were not changed.

Added estimator-specific influence covariograms using the location row of the
full joint bread, including nuisance-scale effects. The scalar influence kernel
quadratic form was verified against the joint HAC location variance. Missing
initial positive covariance or an absent zero crossing remain explicit rather
than receiving a default cutoff. The exponential decay parameter is not a
population zero-crossing distance: that covariance remains positive at finite
lags. No current interval uses an automatically selected influence bandwidth.

A translation-invariance test caught floating-point lag-bin edge assignment.
Corrected right-closed bins and lag cutoffs using coordinate-precision tolerances;
verified identical pair counts/covariances after shifts up to 1e6. Re-ran the
study from identical seeds in a new bins-checked folder, preserving the original.
Critical values, all interval rows, Gaussian validation, and censored coverage
summaries were byte-identical before/after the diagnostic repair. Only the
influence-range diagnostics changed. Use the bins-checked report below.

- [Paper fetch list](PAPER_FETCH_LIST.md)
- [Corrected calibration report](output/gaussian_hac_calibration_bins_checked_2026-10-31/README.md)
- [Censored confirmation summaries](output/gaussian_hac_calibration_bins_checked_2026-10-31/censored_summary.csv)
- [Influence-range diagnostics](output/gaussian_hac_calibration_bins_checked_2026-10-31/influence_range_summary.csv)
- [Study runner](scripts/analysis/validate_gaussian_hac_calibration.R)
- [Influence identity/translation tests](scripts/analysis/verify_influence_covariogram.R)
- [Report and reconciliation](scripts/analysis/report_gaussian_hac_calibration.R)

Saved a source snapshot and SHA-256 list with the corrected study. Tightened the
runner's checkpoint guard to hash study/helper and package R sources and record
R/platform/survival versions. Legacy or changed-source checkpoints are rejected
with a fresh-output-directory instruction; neither caches nor completed results
are deleted. The completed study's snapshot records the checked source before
this final guard tightening; future reruns must use a new directory. Integration
checks confirmed rejection of both legacy and changed-source cache reuse.

Next scientifically justified work: read the fixed-b and spatial-robust inference
papers, assess a theory-supported calibration and bandwidth rule with independent
validation, and distinguish observed-domain uncertainty from informative omitted-
site bias. A favorable conditional simulation is not a universal calibrated
confidence interval. No production API, full-data estimate, dependency, version,
release, or external communication was changed by this research phase.

## Supplied-paper methods review — November 1, 2026

Reviewed the original nine supplied papers and three additional papers/versions
that arrived during reading. The inventory at this milestone contains fourteen
PDF files representing twelve distinct papers/versions, plus second copies of
Ordoñez et al. (2018) and Qiang & Craigmile (2026). Full extracted-text comparisons
confirmed the extra copies (ignoring publisher download stamps for Qiang).
Original PDFs and filenames were preserved.

Saved PAPER_REVIEW.md with paper-specific page references, assumptions,
applicability to the current censored/declustered estimator, limitations,
supplement requests, and the proposed next validation sequence. This was a
methods-focused first pass, not a complete proof audit. Source fingerprints,
PDF page counts, pages consulted, review status, and duplicate comparison are
in output/literature_review_2026-11-01/. Selected method pages were also rendered
to check equations visually. PAPER_FETCH_LIST.md now distinguishes received,
reviewed, earlier-version, and outstanding material.

The review reinforces separating spatial covariance calibration, exogenous
nonuniform sampling, and response-dependent sampling. Conley/fixed-b assumptions
do not automatically cover preferential locations or estimated declustering
weights. Müller–Watson SCPC is a distinct published variance-and-critical-value
procedure: our earlier Gaussian-reference Bartlett study is not SCPC. Lahiri–Zhu
GBBB requires more than square blocks, including its source-block construction
and centered bootstrap equations. The current occupied-block bootstrap remains
an experimental sensitivity method. Diggle's discussion also cautions against
calling every adaptive sampling history nonignorable without considering the
conditioning information retained.

The next recommended experiment is an independent small-design reproduction of
published SCPC on uncensored Gaussian observations, checked against author
material, before trying fitted censored score influences and then estimated
weights. Keep a separate GBBB reproduction and preferential-sampling experiment;
do not tune on earlier confirmation data. No new simulation or production spatial
interval was introduced during this review.

Confirmed exact supplement filenames on publisher pages: Müller & Watson's
ecta200510-sup-0001-dataandprograms.zip and Qiang & Craigmile's
env70081-sup-0001-SupplementalMaterial.pdf (Data S1). The extra Qiang main-paper
copy is not that supplement. Bester (2016) appendices, Bester (2011) nonlinear
appendix, and Diggle discussion/rejoinder are already in the supplied PDFs.
Kelly WP19/13 is the 2019 sole-author working paper, not the requested Conley &
Kelly (2025) published version. The supplied sandwich-OOP document identifies
itself as a modified version of Zeileis's article. These versions are recorded
explicitly rather than silently substituted.

Added build exclusions for Papers/ and PAPER_REVIEW.md. No package R code,
dependency, version, full-data estimate, release, or external communication was
changed. Packaging verification is recorded alongside this review; full package
checks and simulations are not claimed to have been repeated for a literature
and build-exclusion update.

A staged source build with --no-build-vignettes --no-manual passed. All current
reference PDFs and review/fetch notes were present in staging; the resulting
archive contained none of Papers/, PAPER_REVIEW.md, or PAPER_FETCH_LIST.md.
The direct whole-workspace build was interrupted while staging; temporary partial
files were preserved. Final source fingerprints, local review links, duplicate
text comparisons, and git diff whitespace checks passed.

## Supplemental review and Gaussian projection checks — November 1, 2026

Received and reviewed Qiang & Craigmile's eight-page Data S1 PDF and inspected
Müller & Watson's replication ZIP. Updated PAPER_REVIEW.md and
PAPER_FETCH_LIST.md so these supplements no longer appear outstanding. Saved
source hashes, the archive's 17 file-member hashes, review scope, findings,
and limitations in output/supplement_review_2026-11-01/.

Qiang's supplement makes the range-prior assumptions explicit: simulation prior
means match the generating spatial ranges. Its parameter-recovery tables are
individual examples, not repeated credible-interval coverage experiments.
Table S5's pooled LGCP preferential parameter interval includes zero. Future
preferential-sampling validation should vary prior centering/strength and
evaluate repeated coverage separately from prediction. All eight pages were
read; Table S5 was also visually verified after rendering.

The Müller–Watson archive contains Fortran source and state/light data. Its
README specifies Intel Fortran and IMSL and the source contains machine-specific
Windows paths. The original program was not compiled or executed. Read the
main driver, covariance/projection, probability, and calibration modules; other
members were inventoried and dependency/license markers inspected. No author
implementation was copied into the package and no dependency was installed.

Added original, research-only R helpers for fixed-component Gaussian projection
formulas, adaptive quadratic-form integration, and critical-value inversion:
scripts/analysis/scpc_gaussian_reference.R and
scripts/analysis/verify_scpc_gaussian_reference.R. The archive uses a raw-ratio
statistic, whose critical value is the paper's usual statistic critical value
divided by sqrt(q). Verified this convention and the equivalent symmetric
quadratic-form eigenvalues. The source calibration grid has 30 decay values
from c0 to exp(5)*c0; its approximation should not be confused with an exact
continuum supremum.

Verification passed iid Student-t critical values for q=1,2,3,5,10 (maximum
absolute discrepancy below 2e-10), coordinate translation/rotation, projection
orthonormality, and statistic/standard-error identities. Six prespecified
Gaussian checks used 25 locations, fixed q=3, two exogenous spatial designs,
and three specified covariance cases. Critical values were calibrated by
integration, then checked against 100,000 fresh full-vector Gaussian draws per
case: 600,000 total. Rejection frequencies ranged from 4.922% to 5.131%, within
the stated Monte Carlo tolerance of 5%. Results, source hashes, and runtime
information are in output/scpc_gaussian_formula_checks_2026-11-01/.

These checks validate a bounded numerical component, not a complete SCPC
implementation, a rerun of the original Fortran, a covariance-envelope guarantee,
or censored/declustered USGS interval coverage. Component-count selection and
envelope calibration remain the next step, before censored-score adaptation.
No production API, package dependency, version, full-data estimate, or routine
confidence interval changed. Full package checks were not repeated for these
excluded research scripts and notes; the mathematical checks above were run.

## Numerical SCPC envelope and censored transfer — November 1, 2026

Extended the original fixed-component Gaussian formula reference to numerical
covariance-envelope calibration and component selection. New research files:
scpc_envelope_reference.R, validate_scpc_envelope.R,
validate_scpc_censored_transfer.R, and report_scpc_validation.R, all under
scripts/analysis/. No external implementation code was copied.

The Gaussian study uses four fixed designs (36/64 jittered or exogenously
clustered locations), and c0 chosen from locations alone to imply 3% average
pairwise exponential correlation. Candidate q values 1–20 are evaluated using
expected iid interval length. Calibrations use 31 then 61 log-decay nodes with
local peak searches; a separate 241-midpoint grid and local searches audit the
selected statistic. The near-independence tail has an explicit Gaussian
total-variation bound. This is not a global certificate for the compact decay
interval, nor a guarantee under arbitrary covariance or preferential locations.

Selected q/critical values were 10/2.7073 (jittered 36), 7/2.5534 (clustered 36),
8/2.7391 (jittered 64), and 7/2.4752 (clustered 64). Choices were stable under
grid refinement, none reached the candidate boundary, and scale/rotation checks
passed. The largest separately audited rejection probabilities were below 5%.
Calibrations and spatial bases were saved before response simulation.

Independent simulation used 20,000 development plus 100,000 confirmation draws
for each design/decay combination, for 2.4 million full-vector Gaussian draws
and 7.2 million interval evaluations across three comparators. Confirmation
coverage for the numerical-envelope method was 94.990–97.857% within the tested
exponential envelope, and 85.140–88.684% when dependence was deliberately twice
as persistent as allowed. Analytic rejection and simulation agreed within the
prespecified Monte Carlo tolerances. Source snapshots, inputs, component
candidates, calibration summaries, audits, plots, and report are in
output/scpc_envelope_validation_2026-11-01/.

Then tested transfer to estimated censored lognormal location influences with
equal weights and four measurements per site. The latent log-field SD is 0.55,
independent observation noise SD 0.2, left threshold exp(log(50)-0.3), and upper
threshold exp(log(50)+1). The target is the common superpopulation location
log(50), not a realized finite-domain average. Critical values, component counts,
and bases remained frozen from the Gaussian study. All observed sites remain in
the design; no preferential omission or estimated declustering weights were used.

Across four designs and three latent covariance cases, 200 development plus
500 confirmation datasets give 8,400 attempted fits. All fits were available.
Confirmation SCPC transfer coverage ranged 94.6–99.0%; at the envelope boundary
it ranged 94.6–95.6%. The working-independence comparator undercovered even with
independent sites because measurements at a shared site remain correlated.

Added a site-cluster sandwich/t comparator after inspecting that initial
two-method study. This is a secondary comparison, explicitly recorded as such;
it was not used to retune SCPC. The site correction permits within-site
dependence but treats distinct sites as independent. At the envelope boundary,
its confirmation coverage was 77.4–86.4%, compared with SCPC's 94.6–95.6%.
With independent sites it gave 94.2–97.2%. This separates the benefit of handling
site repeats from between-site spatial dependence.

Preserved the initial study in output/scpc_censored_transfer_2026-11-01/ and
reran from identical seeds in a new directory to add the site comparator:
output/scpc_censored_transfer_site_control_2026-11-01/. All 8,400 original result
rows/columns were identical after adding the comparator. The new study has
25,200 interval evaluations, 72 summary rows, source/input fingerprints, snapshots,
failure counts, censoring counts, runtime information, plots, and report.

Reports reconcile all counts; scientific plots were visually inspected. Both
new study runners reject nonempty output directories, and their completed
source/input hashes were checked. Numerical Gaussian checks, independent
censored confirmation, expected iid length, and git diff whitespace checks
passed. Package checks were not repeated: package R code, dependencies, version,
USGS estimates, and routine intervals did not change.

Next step: weighted Gaussian estimating-equation calibration before testing
coordinate-derived declustering weights in censored models, with new independent
designs and confirmation seeds. Preferential sampling and missing support remain
separate modeling problems. The favorable equal-weight result does not validate
a default interval for the historical USGS compilation.

## Conditional declustering-weight calibration — November 1, 2026

Completed the next research stage: weighted Gaussian calibration followed by
independent censored-lognormal transfer. Added scpc_weighted_reference.R,
verify_scpc_weighted_reference.R,
verify_scpc_weighted_censored_reduction.R,
validate_scpc_weighted_gaussian.R, validate_scpc_weighted_censored.R, and
report_scpc_weighted.R under scripts/analysis/. No package R function changed.

The extension explicitly accounts for the weighted estimate and fitted residual
operator. With normalized weights a, theta_hat=a'Y and fitted influences are
psi=B*Y for B=n*diag(a)*(I-1*a'). Components are chosen from B*Sigma(c0)*B';
calibration uses the joint Gaussian transform [sqrt(n)*a, B'*V]. Even at
independence the resulting covariance is generally not identity, so ordinary
Student-t quantiles and the equal-weight expected-length formula are not reused.
Expected iid length is evaluated by a Gaussian-norm integral. This is an
original conditional-Gaussian research extension, not a cited theorem granting
weighted or preferential-sampling validity.

Verification passed weighted residual/mean algebra, equal-weight reduction,
common weight-scale invariance, chi-norm identities and independent norm
simulation. A separate uncensored survreg check reproduced the exact weighted
mean and fitted influence contributions, including invariance to a concentration
unit multiplier and common weight rescaling. The full joint location/log-scale
bread remains in the censored influence calculation.

Four fresh exogenous spatial patterns have 36/64 jittered or clustered sites.
Cell width 0.25 in unit-domain coordinates and 100 origins were fixed in advance.
Weights calculated for four colocated records per site aggregate to the same
site weights. Kish weight n was 35.68/60.68 for the jittered designs and
17.16/24.44 for the clustered designs; this measures imbalance, not independent
replication. Locations and coordinate-derived weights remain fixed conditional
on each design throughout the experiment.

Searched all feasible nonzero component counts, 1 through n-1, with the stated
expected iid interval-length criterion. Selected counts were 12/24/12/33 for
jittered36/clustered36/jittered64/clustered64. The baseline c0 implies 3% average
pairwise exponential correlation. Counts and critical values were stable between
31- and 61-node calibration grids with local searches. Separate 241-midpoint
audits and local searches gave rejection maxima below 5%; the independence tail
has the Gaussian total-variation bound. The compact interval remains a numerical
search, not a global continuum certificate.

The Gaussian study used 20,000 development plus 100,000 independently seeded
confirmation vectors per four-design/four-decay combination: 1.92 million
full-vector draws and 5.76 million interval evaluations across three methods.
Weight-aware coverage within the tested envelope was 94.983–98.890%; under
stronger-than-allowed dependence it was 83.622–92.747%. Analytic probabilities
and Monte Carlo results agreed within the prespecified tolerance; independent
expected interval length also agreed. At the envelope boundary, the two
clustered cases gave 95.125%/94.999% coverage with weight-aware calibration,
versus 67.935%/61.171% when spatial calibration ignored the weights. Both
comparators still use the same weighted point estimate and weighted influences.

Frozen calibrations were then tested on 8,400 new censored fits (200 development
and 500 confirmation datasets per four-design/three-decay combination). Each
site has four measurements with latent log SD 0.55 and observation-noise SD 0.2;
left/upper thresholds remain exp(log(50)-0.3)/exp(log(50)+1). The common target
is the superpopulation log location log(50). No weight, component, or covariance
class was selected or retuned from response data. All 8,400 fits were available.

Weight-aware censored confirmation coverage ranged 95.8–99.6%, with conservative
cases explicitly retained. In the clustered boundary cases it was 96.2%/95.8%,
versus 68.6%/64.0% for geometry-only calibration. Across all clustered cases the
geometry-only comparator ranged 42.0–68.6%, including failure at independent
sites because the weighting transformation was ignored. The site-cluster-t
comparison allows site repeats and weight imbalance but treats distinct sites
as independent; its results are retained separately. There are 25,200 censored
interval evaluations and 72 summary rows.

Reports, same-axis scientific plots, raw results, source snapshots, runtime
information, calibrated weights/bases, hashes, audits, and summaries are in:

- output/scpc_weighted_gaussian_2026-11-01/
- output/scpc_weighted_censored_2026-11-01/

The weighted linear-Gaussian transformation is exact; its transfer to nonlinear
censored influences is empirical. Four fixed location patterns do not establish
unconditional coverage when locations and weights change across surveys.
Preferential sampling, trends, support/method mixtures, unsampled terrain, and
more severe censoring are not validated. The next experiment should recalculate
weights across independently generated site patterns and assess those remaining
assumptions separately. Routine package intervals, full-data USGS estimates,
dependencies, and version remain unchanged; no release or external message was
made. Full package checks were not repeated for these excluded research files.

Final verification reconciled Gaussian/censored counts, source and calibration
fingerprints, report links, and whitespace checks. Both weighted study runners
correctly rejected reuse of their completed output directories. Scientific plots
were rendered and inspected; coverage axes are shared across the four panels.

## Changing survey layouts: frozen protocol — November 1, 2026

Started an independent random-layout extension of the weighted study. The
prespecified design has 20 development and 80 confirmation layouts per two
exogenous sampling regimes, for 200 independent layouts with 24 sites each.
Uniform locations cover the unit square. The clustered design draws 18 sites
from [0,0.25] x [0,0.25] and six from the whole square, retaining positive
sampling support throughout the window. Locations do not depend on chemistry.

Every layout recalculates cell weights (width 0.25, 100 origins), all feasible
component counts 1–23, and calibration/audits before responses are generated.
The physical benchmark decay is fixed at c0=10 across surveys; actual decays
are also fixed physical values, not reset to enforce a geometry-specific
average correlation. Conditional Gaussian validation uses decays 5 (outside
the envelope), 10, 20, and independence, with 1,000 fields per layout/case.

For censored transfer, 25 fresh latent fields per layout/decay (10,20,iid) each
give four observations per site. Each field is paired across moderate and
heavy censoring: left offsets -0.3 and +0.6 from log(50), upper offset +1,
latent SD0.55, record noise SD0.2. This gives 15,000 distinct underlying fields
and 30,000 attempted censored analyses. The target remains the common
superpopulation log location, not a realized finite-domain spatial average.

No exact observations with these common censoring bounds is prespecified as
unavailable: there is no finite two-parameter likelihood maximum in that case.
Other numerical or calibration failures are retained as unavailable intervals.
All attempts count in unconditional coverage; available-only coverage is
reported separately. Summary uncertainty is calculated across independent
layouts, not by counting repeated fields at one geometry as independent surveys.

The original probability routine was vectorized without changing its formula.
Against the prior checked source snapshot, 120 probabilities and 30 critical
values differed by at most 1.67e-16 and 2.49e-14. A separate proposed peak-search
optimization gave identical tested critical values but was slower, so it was
rejected and preserved under optimization_trial/. The original peak search,
statistical method, and study specification remain in use. The interrupted
initial checkpoints resumed with the identical source/runtime/configuration
manifest after restoring that search. No confirmation response was used for
tuning. Sources, per-layout geometry/calibrations, and checkpoints are saved in
output/scpc_random_designs_2026-11-01/.

The study and layout-based report code are validate_scpc_random_designs.R and
report_scpc_random_designs.R under scripts/analysis/. Results below will be
added only after independent confirmation and count/provenance reconciliation.

While the random-layout study ran, identified and reviewed the newly supplied
published CensSpatial software article: Ordoñez, Galarza & Lachos (2024),
SoftwareX 27, 101762, DOI 10.1016/j.softx.2024.101762. All four pages were read;
page 2's function table and SAEM output were visually checked. It is a separate
published reference from the earlier 2021 preprint. Its historical software
metadata says version 3.6; this is not asserted to be today's current version.
The illustration uses 25% left censoring and prediction, not validation of our
heavy-censoring common-parameter interval. Its Supplementary Material, including
section S.2, was added to the fetch list; the exact filename was not confirmed.

The new uWL7AG-Download.pdf is byte-identical to the prior Lahiri–Zhu PDF.
tstat.pdf is another published Ibragimov–Müller copy with matching DOI and no
publisher cover. All files were preserved. Hashes and review scope are in
output/literature_additions_2026-11-01/. CensSpatial and tmvtnorm are not installed
in the checked R library; mvtnorm is available. No library changes were made.

## Changing survey layouts: confirmation results — November 1, 2026

Completed all 200 independent layouts, recalculating declustering weights and
spatial calibration before response generation. All 200 numerical calibrations
were available. The Gaussian check used 800,000 distinct full-vector draws and
1.6 million interval evaluations across the two methods. Confirmation coverage
for the weight-aware interval was 95.030–98.470% within the assumed dependence
envelope. With stronger dependence outside that envelope, coverage fell to
84.414–88.976%; the covariance-class assumption remains consequential.

The censored experiment used 15,000 distinct underlying fields, each paired
across moderate and heavy censoring, for 30,000 attempted fits and 60,000
interval evaluations. Confirmation used 80 independent layouts per sampling
regime and 25 fields per layout/case. Results for the weight-aware interval:

| Confirmation setting | Coverage over all attempts | Available-only coverage |
| --- | --- | --- |
| Moderate censoring, across six regime/decay cases | 96.25–99.05% | 96.25–99.05% |
| Heavy censoring, across six regime/decay cases | 85.50–93.75% | 89.39–93.94% |

Every heavy-censoring case had a layout-based Monte Carlo 95% interval whose
upper endpoint was below 95%. Heavy-censoring mean log-location bias ranged
from -0.142 to -0.082 across confirmation cases, versus -0.010 to -0.003 under
moderate censoring. Failed fits alone do not explain the undercoverage: it
persists among available intervals. This is an identified limitation of the
nonlinear censored transfer, even though its Gaussian benchmark worked within
the stated covariance class. No critical value or estimator was retuned using
these confirmation responses.

Across development and confirmation, 29,817 of 30,000 fits were available.
The 183 unavailable fits comprised 182 datasets with no exact observations
at the two common censoring bounds and one failure to solve the influence
equations. All 15,000 moderate-censoring fits were available. Every attempted
dataset remains in unconditional coverage; available-only results are separate.
This research harness permits a finite fit with one exact observation, while
the public package requires at least two; its availability figures must not be
represented as the public API's behavior.

Uncertainty in measured coverage is estimated across independent layouts.
Fields at a single geometry are not counted as independent surveys. Paired
method and censoring-severity comparisons are retained. Numerical calibration
uses grid refinement, local searches and an independent dense audit, not a
global compact-interval proof. The target is the common superpopulation log
location; this does not validate finite-domain means, preferential sampling,
or actual NGDB confidence intervals.

Report and evidence:

- [Study report](output/scpc_random_designs_2026-11-01/README.md)
- [Confirmation summaries](output/scpc_random_designs_2026-11-01/censored_summary.csv)
- [Paired censoring comparison](output/scpc_random_designs_2026-11-01/paired_censoring_difference.csv)
- [Recorded failures](output/scpc_random_designs_2026-11-01/fit_failures.csv)
- [Final verification](output/scpc_random_designs_2026-11-01/FINAL_VERIFICATION.txt)
- [Artifact fingerprints](output/scpc_random_designs_2026-11-01/final_artifact_manifest.csv)

Final checks reconciled raw and summary counts, source fingerprints and their
snapshots. Resuming the completed study with identical provenance left result
tables byte-identical; a changed-configuration trial was rejected before
checkpoint creation. Both scientific figures were visually inspected. Report
generation and whitespace checks passed. Reporting and probability-verification
source snapshots accompany the original frozen study sources. Full package
checks were not repeated for this continuation's excluded research code and
documentation changes; production R code was unchanged during this study.

Updated package guidance and literature notes to emphasize the heavy-censoring
failure. Added the local Papers directory to Git ignore rules; the supplied
papers remain intact and bibliographic notes remain available. Routine weighted
confidence limits remain withheld. Package version remains 0.1.0 toward the
intended 1.1 milestone; no release, commit, push, or external message was made.

Next priority: compare a spatial-likelihood or bias-aware estimator under heavy
censoring, using a new declared protocol and fresh validation seeds. Assess
point bias and interval coverage separately. The completed confirmation study
must remain evidence of this limitation, not become a tuning set presented
as new confirmation. Preferential sampling requires a separate experiment.

## CensSpatial supplement received and reviewed — November 1, 2026

Ross reports no supplements or code are available for Ordoñez et al. (2018).
Closed that fetch request on his report; the already reviewed main paper remains
available. Found the newly supplied 2024 CensSpatial supplement as
Papers/1-s2.0-S235271102400133X-mmc1.pdf. Reviewed all nine pages and rendered
pages 4, 5 and 7 to check SAEM updates, prediction equations and Figure S.2.1.

The supplement provides an observed censored-Gaussian likelihood, conditional
truncated-normal SAEM construction, and prediction recipes. S.2 is exploratory
rainfall trend/variogram analysis, not a heavy-censoring interval experiment.
Printed equation S.1.9 uses the wrong cross-covariance orientation for the block
dimensions stated above it; this is a dimensional check of the document, not
evidence of a software bug. A model-based likelihood comparator remains the
next proposed benchmark; the supplement does not establish calibrated intervals
for our common parameter or resolve the completed confirmation failure.

Review, page-specific findings, proposed benchmark checks and file fingerprints:
output/censspatial_supplement_review_2026-11-01/. Supplied originals are intact.
The fetch list and paper-review status are updated. No estimator, production R
code, package version or library dependency was changed.

Also identified Lehner's 34-page bandwidth preprint, the four-page geoR newsletter
extract, and a 13-page Fouedjio censored-random-forest paper from their first pages;
full reviews are pending. An 18-page candidate online appendix was inspected on
its first three pages but lacks an internally confirmed title/author identity;
it is not yet counted as the requested Conley–Kelly published article or code.
These limited scopes are explicit in the new source inventory.

## Additional papers: methods review and counterchecks — November 1, 2026

Reviewed the newly supplied Lehner bandwidth preprint, Fouedjio censored-random-
forest article, Ribeiro–Diggle geoR newsletter article, and Conley–Kelly online
appendix. Extracted 69 pages across four PDFs; the review records the methods,
tables and figures inspected and does not claim a full proof/plot audit.
Six selected pages were rendered to verify formulas and results. Originals were
preserved and SHA-256 fingerprints recorded. Full notes and review scopes:
output/additional_methods_review_2026-11-01/README.md and source_inventory.csv.

Lehner's method motivates a covariogram diagnostic, but its county-centroid
table reports up to 12.8% rejection for nominal 5% Epanechnikov tests (14.9%
for Bartlett). The printed absolute-tolerance selector at eta=0 differs from
the described sign crossing; fixed bins and a fixed smooth taper need more
care than a claim of consistent estimation at the true support cutoff alone.
An original deterministic 2D spherical-covariance example confirms persistent
attenuation: at support/bandwidth 1.2 on an expanding unit grid, Bartlett and
Epanechnikov limiting expected variance ratios are 0.886667 and 0.905556.
Direct positive-definite covariance matrices match exact finite-grid pair counts.
This is a check of a generic method claim, not a replication of author code or
a complete audit of every theorem assumption.

Fouedjio's paper is closely relevant to geochemical mapping: the scandium example
has 568 locations, about 6% censoring and 26 domain-wide predictors. Its 147
test observations are selected near censored sites, and error comparisons do
not establish prediction-interval or common-parameter confidence coverage.
Random quadratic projection is also not generally conditional Gaussian sampling:
in the scalar positive-normal example, projection puts half the probability
exactly at zero, while truncation has no boundary atom. Independent integrals
confirm means 0.398942 versus 0.797885. This distinguishes feasible realizations
from a calibrated probability distribution without denying the paper's reported
predictive gains. Repeated-site measurement error remains a separate issue.

The geoR article is useful historical groundwork for an uncensored Gaussian
reference, likelihood profiles and parameter uncertainty. It does not establish
censoring or preferential-sampling support. The 18-page online appendix is now
associated with Conley–Kelly through publisher PII/title/authors/DOI metadata and
matching Appendix C in their author-hosted 2024 working paper. Its first three
pages were methods-reviewed and the remaining plot structure inventoried, not
numerically reproduced. The appendix supplies basis/placebo/cluster diagnostics;
the published main PDF and code are still not present locally.

Added and ran scripts/analysis/verify_literature_method_cautions.R. Deterministic
integration, covariance positive-definiteness, and independent matrix/pair-count
checks passed. These are mathematical checks, not new Monte Carlo coverage
results. Updated the paper review and fetch statuses. No library installation,
production estimator, routine interval, dependency or version change was made;
full package checks were not repeated for excluded research/Markdown changes.
The next benchmark remains the small spatial censored likelihood, with fresh
validation seeds. Mapping and spatial-null regression diagnostics are later
comparisons, not substitutes for addressing the observed heavy-censoring bias.

## Focused spatial likelihood reference: first numerical gate — November 1, 2026

Returned to the measured problem: downward point bias and undercoverage in
heavily censored, clustered data. Implemented an original research-only joint
Gaussian likelihood on log concentration. It partitions exact and censored
responses, conditions the latter on the exact observations, and represents
repeated-site noise explicitly. Covariance is supplied, not estimated. Fixed
integration uniforms are reused across candidate locations and generated
without changing the caller's simulation RNG state. There is no new interval.

The initial reference passed 15 numerical comparisons: complete Gaussian density
and GLS location; independent-censoring reduction; direct conditional integrals;
low-dimensional Miwa comparison; permutation/reflection; a remote univariate
interval; and fixed-uniform reproducibility. RNG preservation, invalid inputs,
singular covariance and no-finite-location cases were checked separately.
Tests and source are excluded research files. Evidence is in
output/spatial_likelihood_reference_2026-11-01/verification.csv.

Declared a small development pilot before response generation: eight new
exogenous layouts, 24 sites and four records each; 16 unique fields paired
across moderate/heavy censoring; 32 attempted analyses. It supplies the true
covariance and compares known-marginal-scale working-independence fits with
equal/cell weights, plus the previous unknown-scale weighted fit. That last
comparison changes nuisance-parameter knowledge as well as spatial modeling.
Fresh seeds, a source/runtime/protocol manifest, checkpoints and independent
integration audits are retained. This is a feasibility pilot, not confirmation.

While it ran, a heavy-censoring numerical failure prompted a separate tail
check. For P(Z1<=10,Z2>=9) under independent standard normals, direct unreflected
CDF subtraction plus a tiny positive integration tolerance produced a finite
but incorrect log probability near -716.714 rather than -43.628. Zero tolerance
exposed the underflow; reflecting the right constraint to the left tail recovered
the analytic value. Numerical details and planned development rerun are in
output/spatial_likelihood_reference_2026-11-01/NUMERICAL_FINDING.md. The original
pilot sources/results stay frozen; no statistical claim is made from finite
optimizer output alone. Tail handling must pass a second gate before scaling.

## Spatial likelihood tail amendment and identification diagnosis — November 1, 2026

The first oracle pilot completed 32 attempts: 31 finite spatial fits and one
numerical failure. Preserved all original data, sources, audits and artifact
fingerprints. Added right-tail/finite-box reflection and zero integration
tolerance to expose underflow rather than return a finite probability floor.
Eighteen numerical equalities plus RNG/input/failure guards then passed,
including correlated remote probabilities against independent 1D integrals.
Checked a mixed-left/right case with no exact observations separately: the
known-covariance location can be finite, unlike the no-exact unknown-scale rule.

A separately declared numerical-development rerun retained the identical eight
layouts, 16 underlying fields, covariance matrices, observed bounds, integration
seeds and comparator fits. All 32 spatial fits were finite. Largest coarse/fine
location shift was 0.003098; largest eight-batch Monte Carlo SE of log likelihood
was 0.022853. The response reuse is explicit: this is not independent
confirmation. A further numerical precision refinement of the largest-MC-error
case was started, with new integration seeds and no new response generation.

Oracle pilot files:
output/spatial_likelihood_reference_2026-11-01/ and
output/spatial_likelihood_reflected_2026-11-01/.
Completed-run resumption preserved the reflected results byte for byte; changed
sources correctly prevented reuse of the original output. Both source snapshots
remain available before subsequent package changes.

Independently diagnosed the finite working-independence fit with the largest
estimated scale (selected post hoc, not confirmation). It had 96 records but
only two exact values carrying 0.922% of the declustering weight. Fitted
log-location was -5.1911 and log-concentration SD 8.1649, versus the generating
marginal SD 0.5852. An independent profile reproduces the maximum; a simple
censoring-bin approximation gives location -5.3545 and SD 8.3034. Almost all
weight is in the two tails of the narrow detection window, so the fit can
match those masses using a huge spread and low mean. Successful convergence
does not establish useful identification. This is distinct from the numerical
integrator failure and cannot be repaired solely by changing a standard error.
No extreme result was clipped. Full profile, fractions, source and plot:
output/censored_scale_identification_2026-11-01/.

## Package improvement: censoring-weight support — November 1, 2026

Added exact_fraction, weighted_exact_fraction, weighted_left_fraction and
weighted_right_fraction to declustering_sensitivity(). The weighted fields are
shares of geometry weight on the same usable chemistry/coordinate subset,
not independent-information counts or inclusion probabilities. They sum to one
when usable rows exist and are missing for empty subsets. They remain available
when the concentration fit is insufficient. Multi-element analysis propagates
them automatically; no new threshold, filter or confidence interval is applied.

Tests cover hand-calculated clustered weights in both censoring directions,
exclusion of invalid chemistry/missing coordinates, all-censored/empty inputs,
and propagation across elements/unsupported units. The full local software test
suite passed. Updated function help, README and NEWS. A source archive built
successfully with vignette creation from a stage honoring .Rbuildignore. The
initial direct build could not copy the mixed local workspace; the staging
approach preserves original files and keeps excluded research artifacts out.
The subsequent full R CMD check and all-USGS support audit were started;
their completion/results will be recorded separately.

The all-data audit reads original chemistry and the prior projected metadata,
recalculates weights for all soil/sediment species at 5/10/20 km with 100 origins,
and requires every usable/excluded/censoring count and Kish weight n to match
the prior all-element results. It refits no concentration estimates and writes
only a new output directory. The 5% comparison uses the existing descriptive
review convention and is not a reliability guarantee.

### Numerical and package checks completed — November 1, 2026

Refined the oracle case with the largest integration audit error using 262,144
and 1,048,576 draws from fresh integration seeds. Locations were 4.060962 and
4.060546, versus the earlier 4.062072. Sixteen independent 65,536-draw audit
batches reduced log-likelihood Monte Carlo SE to 0.006404. Paired likelihood
changes at mu +/-0.01 were both negative with Monte Carlo SE about 0.000114.
This supports a stable known-covariance point calculation at the tested
resolution, not covariance estimation, a global numerical guarantee or coverage.
Evidence: output/spatial_likelihood_precision_2026-11-01/.

The updated package source archive and full local R CMD check --no-manual
completed with Status: OK, including the software suite and vignette creation/
rebuilding. Validation files are in
output/package_validation_weighted_support_2026-11-01/. Archive SHA-256 is
43f70f7a21ff633b53b92c21cc826cd6b39885addb9ad0e74b90752c714d7895.
Installed dependencies were checked; remote CRAN/Bioconductor indexes were
unavailable to the check console. No online dependency-index verification or
PDF manual build is claimed. The archive excludes Papers, research scripts,
generated outputs and the sibling PFAS package. Package version remains 0.1.0;
no commit, push or release was made.

### All-USGS censoring-weight audit completed — November 1, 2026

Completed all 80 soil and 81 sediment species. Every one of the 22,431
element/state/cell rows matched the prior usable/excluded counts, censoring
counts and Kish weight n. There are 11,091 soil rows and 11,340 sediment rows.
Source/input fingerprints stayed unchanged; a completed-run resume preserved
all output tables byte for byte. Public-API spot checks passed on Ag/CT soil
and Ag/DC sediment. The fractions are retained for all rows, including missing
fractions for empty subsets; no concentration estimate or interval was refitted.

Using the existing descriptive review conventions (at least ten exact records
and 5% exact raw share), 322 cell/element/state rows still had under 5% of their
declustering weight on exact values. They represent 131 distinct medium/element/
state combinations across 28 species labels. These overlapping scales are not
independent replicates, and the thresholds are not statistical guarantees.
As one example, South Carolina soil Au has 50 exact values out of 932 usable
records (5.365%); their weight share at 10 km is 0.141%. This shows why raw
detection counts alone can conceal the effect of clustering on the weighted fit.

Report, all retained rows, review cases, plots, counts, checkpoints and hashes:
output/usgs_weighted_censoring_support_2026-11-01/. Both-medium scatterplots were
visually inspected, local report links checked, and no original analysis output
or source chemistry was changed. The full package check remains Status: OK.

### Removing supplied scale: selected-case numerical development — November 1, 2026

Started the next focused gate: estimate a common covariance amplitude while
retaining the true spatial correlation shape and nugget/signal ratio. This
removes known marginal scale but remains an oracle-shape diagnostic, not full
covariance estimation. The reference uses delta=mu/multiplier and
kappa=1/multiplier, with covariance S/kappa^2. Ten numerical equalities passed
against complete Gaussian density, GLS location, exact-data variance MLE,
rebuilt conditional likelihood and independent censored Gaussian fitting.
No-exact non-common finite boxes are explicitly outside the implemented scale
classification rather than incorrectly being declared to have no finite MLE.

Declared four post-hoc selected prior-pilot cases, two coarse and two refined
optimizer starts each, retained bounds/failures, and independent integration
audits before fitting them. Selection targets the prior numerical failure,
largest working-independence scale and its moderate partner, and a one-exact
case (private harness only; public minimum remains two). All reused responses
are development data, not independent confirmation. Protocol and sources:
output/spatial_scale_development_2026-11-01/. No routine interval is introduced.

### Estimated-scale development completed — November 1, 2026

All four selected cases completed with finite solutions from all 16 coarse/
refined optimizer starts; no invalid objective evaluation occurred. Ten prior
numerical equalities and guard checks support the parameter transformation.
Sixteen local parameter perturbations were evaluated with eight independent
integration batches each. Starts, convergence messages, refit shifts and
numerical Monte Carlo errors are retained. Source/input fingerprints matched,
and a completed-run resume left the result table unchanged.

For the paired clustered field, moderate censoring gave log location 3.806456
and marginal SD 0.589605; heavy censoring gave 3.291123 and 0.887701. The
generating values were 3.912023 and 0.585235. The extreme weighted-independent
fit in that heavy case had location -5.191102 and SD 8.164929. The spatial
estimated-scale fit avoided that extreme result but still missed the target.
In the one-exact case it instead estimated location 4.390993 and SD 0.061481,
showing downward scale estimation is also possible. That private one-exact
case remains excluded by the public package's minimum-two-exact fitting rule.

These are deliberately selected development cases, not performance estimates
or independent confirmation. Spatial correlation shape and nugget ratio still
use generating values, so this is not a fully estimated covariance solution.
Finite optimization and negative local perturbation checks do not certify
statistical reliability. The next scientific validation must include estimated
scale explicitly, then unknown range/nugget and fresh response/layout seeds;
it must not reuse these cases as confirmation or introduce routine intervals
from a favorable oracle result.

Report and evidence:

- [Estimated-scale selected-case report](output/spatial_scale_development_2026-11-01/README.md)
- [Every optimizer start](output/spatial_scale_development_2026-11-01/optimization_starts.csv)
- [Independent numerical audits](output/spatial_scale_development_2026-11-01/local_integration_audit.csv)
- [All-USGS censoring-weight audit](output/usgs_weighted_censoring_support_2026-11-01/README.md)
- [Checked package archive and build evidence](output/package_validation_weighted_support_2026-11-01/README.md)

The only production change in this continuation is the descriptive censoring-
weight support in declustering/batch outputs, with documentation and tests.
Experimental likelihood functions remain excluded research scripts. Prior USGS
estimates, all supplied papers and original failed numerical runs are preserved.
Version remains 0.1.0 toward 1.1.0; no commit, push, release, new dependency,
external message or automatic data exclusion was made.

## New geochemical database/sample-size papers — November 1, 2026

Found three new PDFs: Zhang, Manheim, Hinde & Grossman (2005), Statistical
characterization of a large geochemical database and effect of sample size;
Olea (2021), Revisiting the declustering of spatial data with preferential
sampling; and Lawley et al. (2022), Data-driven prospectivity modelling of
sediment-hosted Zn–Pb mineral systems and their critical raw materials.
Extracted 53 pages and completed a methods-focused review, with eight selected
methods/table pages rendered. Page scopes and SHA-256 fingerprints are explicit
in output/geochemical_database_papers_2026-11-01/; originals are intact.

Zhang's work is directly related to the USGS source family, but uses a distinct,
more analytically coherent NGS subset (48,544 samples, reduced to 16,511 by
nearest-grid-node selection at 17 km). It is not the same mixed-method NGDB
best-value compilation or validation of our occupied-cell weights. Its detailed
27-element analysis excluded thirteen elements for detection sensitivity or
precision, including Au/Ag. The paper uses half-LOD substitution for graphical
display. Those conventions are not imported into our likelihood or all-element
outputs. Its 100 repeated selections at sizes 50–5,000 demonstrate test-power
effects; small-subset non-rejection does not make a global distribution correct.
Lilliefors correction materially changes the pass rates and does not correct
spatial dependence, censoring, case weights or rounding. Normality tests are not
generally valid only at small n; no n=1,000 cap, p-value gate or thinning-until-
passing rule is introduced.

The main practical implication is to assess censoring-aware distribution shape
within defined domains, using the existing fitted/nonparametric CDF diagnostics
and scientifically important departures. AIC remains relative model comparison.
Geological mixtures in actual data are an additional adequacy concern; they do
not explain away our heavy-censoring simulation failure under a genuinely
Gaussian log-response or replace validation with estimated scale/range/nugget.

Olea's method repeatedly draws one observation per occupied cell, ranks each
draw, takes rankwise medians and restores missing extrema. Its coal/synthetic
examples are useful comparators but provide no censored-likelihood coverage
validation or empty-cell reconstruction. Original deterministic checks show:
(1) the expected resample mean/CDF equals normalized occupied-cell weighting
for a fixed grid; (2) rankwise median recombination need not preserve balanced
cell representation, even under complete enumeration of equally likely pairs.
This is a geometry/bookkeeping check, not replication of the author's code or
a rejection of its reported example performance. Preserve IDs, geometry and
censoring bounds in any future comparison. The primary USGS software record
confirms a Fortran release including decluster; no code was downloaded or run.

Lawley's prospectivity work contributes geological covariates and spatially
blocked predictive validation. Its AUC concerns deposit ranking/search area,
not calibrated concentration intervals; its evidence weights differ from our
geometry weights. Its imputation/outlier handling is not imported for censored
assays, and a new prospectivity workflow is not started. Supplements were not
reviewed or requested as prerequisites for our current work.

Updated PAPER_REVIEW.md and PAPER_FETCH_LIST.md with receipt/review status.
Ran scripts/analysis/verify_stratified_resampling_interpretation.R; exact
enumeration checks passed. No production R function, package version/dependency,
USGS estimate, filtering threshold or interval changed. Full package checks
were not repeated for reading notes and excluded research code; the previously
checked package remains the current production source.


## Empirical CDF bounds with censoring — November 1, 2026

Implemented the practical distribution-shape diagnostic motivated by the recent
geochemical database/sample-size review. New exported functions
`censored_cdf_bounds()` and `plot_censored_cdf_bounds()` show the minimum and
maximum empirical cumulative fractions compatible with the usable records,
left/right limits, and optional observation weights. No censored concentration
is substituted and no distribution is fitted. All-censored subsets are supported;
empty usable subsets return missing bounds. Positive finite thresholds are
explicit, sorted and unique. Evaluation uses sorted endpoints without allocating
an observations-by-threshold matrix. Weight normalization handles very large
finite weights without overflowing their sum.

The bounds describe observed-record empirical distributions, conditional on
weights and censoring statements. They are not population bounds, confidence
bands, an automatic goodness-of-fit rejection rule, or a correction for missing
terrain/preferential sampling. Fully observed records can have coincident bounds
while sampling and spatial uncertainty remain. No ordinary iid normality-test
p-value, new fit eligibility cutoff, or routine spatial interval was added.
The public minimum-exact-observation fitting rules are unchanged.

Tests include hand-calculated weighted bounds, every combination of nine small
admissible censored-value realizations, exact ties, censoring endpoint inclusion,
all-left/all-right/empty data, zero weights, renamed columns, weight scaling and
replication invariance of the bounds. Right-censored values are strictly above
their limits, so upper-bound inclusion uses limit < threshold. Plotting uses
vertical ranges at evaluated thresholds without interpolating between them.
Updated help, exports, README, NEWS and the worked synthetic USGS vignette.

Ran a selected real-data example on all usable South Carolina soil Au and Al
records, with equal weights and 10 km cell weights averaged over 100 grid
origins. Au has 932 usable records, of which 50 are exact (5.365%); those exact
values carry only 0.141% of cell-declustering weight. Al has 691 usable records,
684 exact (98.987%), carrying 99.209% of cell weight. The resulting plot shows
strong unresolved censoring for Au and nearly coincident bounds for Al. All
weighted censoring shares reconcile to the previous all-USGS support audit.
The examples are selected diagnostics, not independent statistical validation;
a state/medium subset can still mix geology, support and analytical methods.

Evidence and preservation:

- [USGS example report and four-panel plot](output/censored_cdf_diagnostics_2026-11-01/README.md)
- [Threshold-level bounds](output/censored_cdf_diagnostics_2026-11-01/empirical_cdf_bounds.csv)
- [Observation/weight support](output/censored_cdf_diagnostics_2026-11-01/observation_support.csv)
- [Package build/check report](output/package_validation_cdf_bounds_2026-11-01/README.md)

The full local test suite passed. A clean 98-file source stage honoring
.Rbuildignore built and passed R CMD check --no-manual with Status: OK,
including tests and vignette creation/rebuilding. All staged files match current
package source. The archive excludes local papers, outputs, research scripts and
the sibling PFAS package. SHA-256:
aa8dfa13b32586cae42ae8ee6c329830be7cb0e775f54d58d5dd6f9c903f948b.
Installed dependencies were checked; remote CRAN/Bioconductor indexes were
unavailable, so no online index verification is claimed. PDF manual generation
was omitted. git diff --check passed. The demo retains source snapshots, input
fingerprints and numerical outputs; all 29 source/input fingerprints remained
unchanged at final verification. Original source data and prior outputs remain
intact. Version remains 0.1.0 toward 1.1.0; no dependency or release changed.

The next scientific step remains fresh validation with estimated spatial scale,
then unknown range/nugget. These descriptive bounds do not resolve the earlier
heavy-censoring bias/coverage findings or certify experimental spatial intervals.


## Fresh estimated-scale pilot initiated — November 1, 2026

Following the empirical CDF diagnostic, froze a new point-estimation protocol
before generating fresh responses. The study uses 16 independently generated
exogenous layouts (eight uniform/eight clustered), 24 sites and two records per
site. Two spatial decays and paired moderate/heavy censoring produce 32 fields
and 64 attempted analyses. The two-record design is explicitly different from
the earlier four-record selected cases; this is a small fresh pilot, not a
coverage or national-scale validation. Spatial shape/nugget ratio remain known;
one common covariance scale is estimated. Comparators retain known-scale
spatial, equal-independent, weighted-independent and fully observed benchmarks.

Protocol: output/spatial_scale_fresh_2026-11-01/PROTOCOL.md.
Driver: scripts/analysis/validate_spatial_estimated_scale_fresh.R.
All source/protocol fingerprints and runtime are frozen in the study manifest.
Unavailable fits and independent numerical-audit flags will be preserved, with
sampling Monte Carlo uncertainty grouped by independent layout. No production
package code, fitting eligibility, existing USGS output or dependency changed.
Results and interpretation will be appended after the frozen run completes.


### Database terminology check — November 1, 2026

Verified the 1.52 million figure against the USGS release description,
Geochemical data for rock, sediment, soil, and mineral samples, United States
and territories, 1962 to 2023: it reports 1,524,014 samples, not a verified count
of unique sites. Source: [USGS data catalog](https://data.usgs.gov/datacatalog/data/USGS%3A63a107a0d34e0de3a1f2795f).
This full-archive count must be distinguished from the soil/sediment subsets
used in our existing local analyses. Lower/upper reporting limit is preferable
to calling upper censoring a detection limit. Current production parsing handles
left/right censoring; it does not yet decode semi-quantitative reporting classes
into finite intervals. Wording that includes that capability must identify it
as planned work, not a completed feature.


## Fresh estimated-scale pilot completed — November 1, 2026

Completed the frozen 16-layout, 32-field, 64-analysis pilot. Every survey has
24 sites and two records/site. Uniform and clustered geometry is exogenous;
moderate/heavy censoring is paired within each field. The spatial shape and
nugget/signal ratio remain known; a common scale is fitted. This is a fresh
small-design point-estimation pilot, not coverage validation or replication of
the previous four-record selected-case design.

All 64 spatial fits and 256 coarse/fine optimizer starts were finite, with
zero invalid probes and zero prespecified numerical-review flags. Maximum
absolute coarse/fine location shift was 0.005677, log multiplier shift 0.006713,
and central loglik integration MCSE 0.008237. Every start and eight-batch local
audit is retained. Numerical stability does not certify statistical reliability.

Moderate-censoring location RMSE was 0.194 (uniform) and 0.163 (clustered).
Heavy-censoring RMSE was 0.646 and 0.536, compared with 0.202 and 0.207 when
true scale was supplied. However, both largest estimated-scale errors occurred
in private one-exact cases, which the public fitting rule already excludes.
A clearly labeled secondary descriptive restriction to the existing two-exact
count rule leaves 62 cases; heavy-censoring RMSE becomes 0.307 and 0.215.
Both primary and secondary summaries are retained. Do not present the primary
extremes as outcomes of public-package-eligible fits or infer that the existing
count rule alone establishes reliability. A two-exact case still estimated
location 2.996 (geometric mean about 20) versus true location log(50).

Signed heavy-censoring biases -0.089/-0.126 have layout MCSE 0.168/0.140;
this small study does not establish systematic downward bias. Observed errors
were smaller than working-independence comparators overall, but common-case
paired differences are imprecise. Sampling MC uncertainty is grouped across
eight independent layouts/regime; numerical integration MC error is separate.
No 64-independent-survey, model-dominance, interval-coverage or fully estimated
covariance claim is made. Unavailable fits were prespecified for retention;
none occurred here.

Evidence:

- [Scientific interpretation](output/spatial_scale_fresh_2026-11-01/INTERPRETATION.md)
- [Protocol, full tables and plots](output/spatial_scale_fresh_2026-11-01/README.md)
- [Every attempted case](output/spatial_scale_fresh_2026-11-01/fresh_results.csv)
- [Primary method summaries](output/spatial_scale_fresh_2026-11-01/performance_summary.csv)
- [Secondary existing-count-rule summaries](output/spatial_scale_fresh_2026-11-01/secondary_count_eligible_summary.csv)
- [Paired method comparisons](output/spatial_scale_fresh_2026-11-01/paired_comparisons.csv)

Driver/protocol/runtime fingerprints matched at completion and on a completed
resume. The completed resume left all 161 then-existing study artifacts byte-
identical. Verified 16 unique layout seeds, 32 unique field seeds, two records
per site and normalized weights; a saved uniform layout reproduced exactly in
an independent R process. Analytic and numerical complete-data benchmarks agree.
Report MCSE checks cover balanced/unequal cluster counts and one available
layout. Both plots were rendered and visually inspected. Initial report source
is preserved alongside the final source; the reporting-only update removes a
deprecated plotting operation and adds the labeled secondary count analysis.
No likelihood, fit, source dataset or frozen outcome changed after inspection.

Only excluded research scripts and development records changed this turn.
The 98-file previously checked package stage still matches production source,
so package checks were not redundantly repeated. Version remains 0.1.0 toward
1.1.0. No new dependency, public API, fitting threshold, interval, USGS estimate,
commit, push or release was introduced. Original papers, prior numerical failures
and all previous analysis outputs are preserved.

Next: inspect the joint location/scale likelihood in count-eligible weak cases,
explicitly treating their selection as post hoc development. This should guide
an identification diagnostic before unknown range/nugget fitting. Any resulting
method change needs a new prespecified validation with fresh responses; this
pilot cannot serve as its confirmation set.


## Selected spatial scale profiles initiated — November 1, 2026

Following the fresh estimated-scale pilot, froze a post hoc mechanism diagnostic
on the largest count-eligible heavy-censoring location error in each geometry
regime, plus its moderate partner: uniform_06_decay10 and
clustered_07_decay20. These are two existing fields/four censored cases, not
independent confirmation or a sample of failure frequency. Covariance shape
and nugget/signal ratio remain supplied at their generating values.

New excluded research helper profiles covariance scale at a fixed location,
records all scalar evaluations/unavailable probes, and labels search-boundary
or grid-only outcomes instead of treating numerical bounds as estimates.
Before calculating cases, ten analytic/deterministic comparisons passed:
fixed-location complete Gaussian scale MLE/density at four locations and a
censored simulation-vs-Miwa example. Boundary, argument and RNG guards passed.
The helper does not add a confidence cutoff or public eligibility threshold.

Protocol: output/spatial_scale_profiles_2026-11-01/PROTOCOL.md.
Profiles use common fixed uniforms, new-uniform anchor refinements and eight
independent batches evaluating frozen anchor parameters. No re-optimization
occurs within numerical audit batches. Results will distinguish numerical
Monte Carlo error from finite-sample likelihood support; no interval coverage
claim follows. Driver, helper, inputs and protocol fingerprints are frozen.


## Selected spatial scale profiles completed — November 1, 2026

Resumed the saved diagnostic and found all four cases had finished. Generated
and visually inspected likelihood and profiled-scale plots. This is post hoc
mechanism development on two existing fields: the largest count-eligible heavy-
censoring error in each geometry regime plus its moderate partner. There are
not four independent surveys, and no failure-frequency or coverage inference
follows from this selection. Both heavy cases already trigger the existing
default ten-exact review convention despite meeting the two-exact fitting rule.

The uniform heavy case has 45 left-censored, two exact and one right-censored
record, with the exact observations at two distinct sites. Its refitted location
is 2.996617 (geometric mean about 20.02), with marginal log-response SD 1.072828.
Moving location down by one, to a geometric mean about 7.36, costs only 0.222965
loglik units after refitting scale; the paired numerical MCSE is 0.000214.
The upward one-unit loss is 1.507693. This is broad asymmetric finite-sample
support, not evidence that either mean has a particular posterior probability.
The true generating geometric mean is 50, but the examples were selected using
that truth and cannot become a method-confirmation set.

At one log unit below each fitted location, heavy/moderate losses are
0.223/12.854 for the uniform field and 0.669/5.429 for the clustered field.
These offsets are relative to each fit, not the same absolute location across
conditions. The clustered heavy case has four exact records at four distinct
sites, so duplicate exact records alone do not explain these examples.
The fixed-scale and profiled curves share the free-model reference maximum;
fixing generating scale supplies extra knowledge and is not a practical remedy.

All 324 main location profiles have finite interior scale solutions, with no
search-boundary or grid-only outcomes. Deliberately extreme trial parameters
produced 1,231 unavailable likelihood probes, all preserved explicitly; these
are not failed profile maxima and were not clipped to finite probability floors.
All eight joint-reference starts are finite with zero invalid probes. Maximum
source/refit location shift is 0.000419; maximum absolute log multiplier change
across 24 new-uniform anchor profiles is 0.000425. Eight independent integration
batches evaluate frozen anchor parameters and retain paired numerical MCSEs.
Ten analytic/deterministic comparisons and boundary/argument/RNG guards passed
before the selected calculations. Numerical errors are distinct from sampling
uncertainty and do not certify a global optimum or population interval.

Evidence:

- [Scientific interpretation](output/spatial_scale_profiles_2026-11-01/INTERPRETATION.md)
- [Full report and figures](output/spatial_scale_profiles_2026-11-01/README.md)
- [All location profiles](output/spatial_scale_profiles_2026-11-01/location_profiles.csv)
- [All scalar evaluations/statuses](output/spatial_scale_profiles_2026-11-01/all_profile_evaluations.csv)
- [Independent anchor audits](output/spatial_scale_profiles_2026-11-01/independent_anchor_audit.csv)
- [Exact-site support](output/spatial_scale_profiles_2026-11-01/exact_site_support.csv)
- [Archive verification scope](output/spatial_scale_profiles_2026-11-01/ARCHIVE_VERIFICATION_NOTE.md)

All 16 frozen source/input fingerprints match. Completed-driver resume preserved
all 28 then-existing diagnostic artifacts byte-for-byte. The previous fresh-pilot
manifest matches all 85 locally available files; 79 older cloud-only files were
not reread. Optional wider hash comparisons waited on cloud-only archived files
and were interrupted without changing them. The full old staging-tree/current-
source comparison was not completed in this continuation; no fresh full-archive
verification is claimed. New diagnostic inputs, results and figures are local
and verified. Final sources and artifact hashes are retained separately.

Only excluded research scripts and development records changed. No production
package function, confidence interval, count threshold, dependency/version,
USGS estimate, release, commit or push was changed. The prior package check is
still documented; it was not redundantly rerun for excluded research work.
Original papers, source records and previous numerical failures are preserved.

Next gate: evaluate estimated covariance parameters and uncertainty on fresh
prespecified data. Any identification/interval rule informed by these selected
profiles requires new validation. Broad finite-sample support is distinct from
structural non-identifiability, and no Wilks cutoff or automatic filter is
introduced from these examples.

Targeted whitespace checks passed for the four new research scripts and this
log. A repository-wide git diff --check waited on older cloud-only files and
was interrupted; no fresh repository-wide whitespace result is claimed.


## Fresh estimated-covariance pilot initiated — November 1, 2026

User reaffirmed "No data = No data" and authorized fresh validation with
covariance parameters estimated from observations; confidence limits remain
experimental. No concentrations are generated for unsampled terrain, and
censoring limits remain observed bounds rather than substituted values.
Unresolved accepted parameter fields remain missing, with separate raw
candidates and diagnostics retained.

New excluded research code estimates common log location, marginal SD,
exponential range and record-noise/nugget fraction together from observed
bounds and coordinates. It receives no true covariance, generating parameters,
latent responses or oracle control estimates. The highest finite candidate is
selected before convergence/search-boundary classification; a stronger
unresolved candidate is not replaced by a lower-likelihood interior one.
All outside-search and unavailable likelihood probes are recorded without a
finite likelihood floor. No-data and no-exact scale-identification failures
have explicit missing estimates.

Before fresh responses, ten numerical comparisons passed, including covariance
identities, exact multivariate Gaussian/independent censored densities, a
simulation-vs-Miwa rectangle, and agreement of the four-parameter complete-data
fit with a separate two-dimensional optimization that profiles location/scale
analytically. Boundary, unresolved-candidate, no-data/no-exact, singular-matrix
and RNG guards passed. The verification field is development only.

Frozen pilot protocol: output/spatial_covariance_fresh_2026-11-01/PROTOCOL.md.
Eight new exogenous layouts/fields (four uniform/four clustered), 24 sites,
two records/site and paired moderate/heavy censoring give 16 analyses, not
16 independent surveys. Spatial range and nugget ratio are estimated within
a specified exponential-plus-noise family; this is not unrestricted covariance
validation. Controls on the same observations retain a known-shape spatial fit
and an independent fit, separately labeled and never used to initialize the
primary estimator. Source/protocol/runtime hashes and case checkpoints are
frozen. No public API, data filter, production interval or version changed.


## Fresh estimated-covariance pilot and numerical follow-up completed — November 1, 2026

Completed eight new exogenous layouts/fields, four uniform and four clustered,
with 24 sites and two records/site. Paired moderate/heavy limits produce sixteen
analyses on eight independent surveys. The estimator learns common log location,
marginal SD, exponential range and nugget fraction from observed bounds and
coordinates. It receives no generating covariance, true parameters, latent
responses or oracle estimates. The covariance family remains specified;
this is a small numerical/point-estimation pilot, not unrestricted covariance
recovery, preferential-sampling correction or confidence-interval validation.

The frozen primary results contain thirteen finite-interior fits and three
search-boundary outcomes, all boundaries under heavy censoring. One of the
thirteen interior fits has a numerical-review flag. All eight moderate cases
are interior without numerical flags. All sixteen cases meet the existing
minimum-two-exact count rule; every heavy case has fewer than ten exact records
and would receive the existing default few-exact review warning. No count rule
was changed. Missing observations and unsampled terrain were not filled;
censored records retain their observed bounds.

Primary unresolved accepted location/SD/range/nugget values remain NA. Separate
raw candidates preserve the fitted search-limit values without promoting them
as estimates. The three primary boundary candidates approach the lower nugget,
upper nugget and lower range search limits, respectively. The highest finite
candidate controls classification before boundary/convergence acceptance; a
weaker interior candidate is not substituted for a stronger unresolved one.
A search-boundary result is an unresolved numerical fit, not a claim that no
observed data or no statistical MLE exists.

All eighty primary optimizer candidates report convergence, yet the above
boundary and numerical warnings persist. There are 15,096 recorded likelihood
probes: 13,947 finite and 1,149 explicitly outside the search, with no unavailable
within-search probe. Outside trials return Inf as rejected objectives, never a
finite likelihood/probability floor. All sixteen independent candidate audits
are available; maximum central loglik integration MCSE is 0.006949. The ten
pre-response numerical comparisons and missing/boundary/RNG guards passed.
Full optimizer traces, all attempts and all control denominators are retained.

Available-only location RMSE is 0.201 (uniform/moderate), 0.077 (clustered/
moderate), 0.398 (uniform/heavy, only two interior fits) and 0.173 (clustered/
heavy, three initially interior fits including the flagged one). These are
conditional descriptive summaries, not superiority claims; controls have
different available subsets. Common-case comparisons and sampling MCSEs grouped
by four layouts/regime are explicit. Numerical MCSE is separately labeled.
No population confidence interval or likelihood-ratio cutoff was calculated.

### Separate numerical follow-up and endpoint diagnostic

Selected every primary finite-interior case carrying a numerical flag; only
clustered_01_heavy qualifies. This post hoc follow-up uses 32,768 new uniforms,
maxit 2,400 and tighter reltol 1e-8, with two data-only starts. Both integration
fidelity and optimization controls change, so the result is not attributed
solely to Monte Carlo precision. The primary outcomes remain unchanged.

The refined fit is search-boundary: raw location 3.863040, marginal SD 0.377884,
range 0.041585 and nugget fraction 0.979998 near the numerical upper limit 0.98.
Its accepted parameter fields are NA. Location shifts only -0.000857 from the
primary candidate; the covariance classification is the issue. Eight independent
32,768-draw audits retain the improving nugget direction beyond the numerical
search limit (loglik gain 0.0000111, numerical MCSE 0.00000197).

A deliberately separate post hoc independent-noise endpoint check uses the
observed-data independent fit already retained in the primary study. At nugget
fraction one, correlation range drops out exactly and is recorded as NA, not a
placeholder. Endpoint minus bounded-candidate loglik is 0.0001163 with numerical
MCSE 0.0000544: nearly equivalent fitted likelihoods, not a test establishing
independence or a new model-selection rule. The endpoint was outside the primary
finite search box and does not replace its missing parameters. A stable mean
and an unsupported range can coexist; unresolved covariance does not mean all
location information is absent.

The original primary table remains thirteen interior/three boundary. After
follow-up, twelve primary cases are unflagged interior outcomes and four cases
remain unresolved by these numerical searches (three primary boundaries plus
the refined flagged case). These counts do not recompute primary performance
or certify the remaining fits. Final plots explicitly mark the flagged point;
the initial report source/README are preserved alongside the updated report.

Evidence:

- [Fresh covariance interpretation](output/spatial_covariance_fresh_2026-11-01/INTERPRETATION.md)
- [Primary report and figures](output/spatial_covariance_fresh_2026-11-01/README.md)
- [Every original attempt and missing values](output/spatial_covariance_fresh_2026-11-01/covariance_results.csv)
- [Method availability and conditional errors](output/spatial_covariance_fresh_2026-11-01/availability_and_errors.csv)
- [Independent primary numerical audits](output/spatial_covariance_fresh_2026-11-01/independent_numerical_audit.csv)
- [Higher-precision follow-up](output/spatial_covariance_precision_2026-11-01/README.md)
- [Refined comparison](output/spatial_covariance_precision_2026-11-01/refinement_results.csv)
- [Independent endpoint with missing range](output/spatial_covariance_precision_2026-11-01/independent_endpoint.csv)

Primary CSV reconstructs byte-for-byte from all sixteen saved case rows.
Verified paired censoring, two records/site and all eight generating covariance
matrices against a separately assembled parameterization. Completed-driver
resumes preserved all 63 then-existing primary artifacts and all 15 follow-up
artifacts byte-for-byte. Frozen source/input/runtime fingerprints match; final
artifact manifests and report snapshots are retained. Both final figures were
visually inspected. Only excluded research scripts and development records
changed; package checks were not redundantly repeated. No public API, dependency,
version, real USGS estimate, release, commit or push changed. Original records,
papers and prior failed numerical evidence remain intact.

Next development should represent valid covariance endpoints explicitly, with
missing values for parameters that drop out (for example, range at the noise-
only endpoint). Site-only/zero-cross-site-range and noise-only limits need
separate treatment before interpreting finite search caps scientifically. A
revised estimator requires another frozen study with fresh responses; this
pilot and selected follow-up are development data for that revision. Confidence
limits remain experimental.


## Explicit covariance endpoints initiated — November 1, 2026

Started the next research revision after the estimated-covariance pilot: fit
finite positive-range spatial covariance, the site-only/zero-cross-site-range
limit and independent record noise explicitly. Positive spatial range is NA
in either endpoint branch. Candidate likelihoods are compared on independent
integration batches; numerically overlapping comparisons stay unresolved with
missing overall parameters. This is numerical resolution, not family testing,
confidence calibration or automatic fallback to the simplest branch. Stronger
unresolved candidates are not replaced by weaker resolved ones.

Before fresh responses, deterministic checks exposed two development issues.
The no-exact classification was extended to equal reporting limits (sign-only
observations still do not identify scale). Miwa128 was then inaccurate for a
near-perfectly correlated pair: at nugget .0002, the zero-threshold quadrant
log-probability error was -0.0002339. The failed code and findings are preserved
in output/spatial_endpoints_fresh_2026-11-01/development_snapshot and
NUMERICAL_DEVELOPMENT.md. Revised one-sided bivariate probabilities use TVPACK;
final site candidates are independently checked with Miwa4096. No fresh response
was generated before these fixes and eleven analytic/factorization comparisons
plus endpoint/tie/missing/boundary guards passed.

Frozen endpoint pilot: four new exogenous layouts, twelve fields (spatial,
site-only and independent-noise on each layout), twenty-four paired severity
analyses, 24 sites and two records/site. Only two layouts per geometry regime
make this an endpoint-behavior/numerical pilot, not general performance or
coverage evidence. Two independent R processes use frozen case seeds and
separate field checkpoints; initialization hashes every generated field/layout.
No generating label, covariance or latent response is passed to the estimator.
The protocol, sources, numerical failure evidence and runtime are fingerprinted.
No public API, confidence interval, data fill, threshold or release changed.


## Explicit covariance endpoints completed on fresh data — November 1, 2026

Implemented independent record noise and the site-only/zero-cross-site-range
limit explicitly alongside finite positive-range spatial covariance. Endpoints
have no estimated positive-distance range; their range fields are NA. Site-only
likelihood factors by exact-coordinate groups, preserving both records: exact
pair density, conditional univariate mass or bivariate one-sided probability.
There is no near-neighbor threshold, record collapse or missing concentration
fill. The independent endpoint fixes nugget fraction structurally to one while
estimating location/SD from observed bounds. Zero-noise singular and arbitrary
long-range limits remain outside this initial endpoint treatment.

Before any validation response, eleven analytic/factorization comparisons passed,
including high-correlation quadrant probabilities, complete Gaussian block
likelihood, a separate balanced site variance-component MLE and stable tail
intervals. Endpoint/tie/missing/stronger-unresolved guards passed. The initial
no-exact equal-threshold classification issue and inaccurate Miwa128 near-
perfect-correlation calculation were corrected before fresh response generation;
failed source and numerical findings remain in development_snapshot and
NUMERICAL_DEVELOPMENT.md. TVPACK is now used for one-sided bivariate site blocks;
final site candidates are independently checked with Miwa4096. Old studies and
their frozen sources were not changed.

The frozen new study uses four exogenous layouts (two uniform/two clustered),
twelve fields (spatial, site-only and independent noise on each layout), and
paired moderate/heavy censoring: twenty-four analyses, not twenty-four
independent surveys. Each field has 24 sites and two records/site. The estimator
receives only observed bounds, coordinates and integration seeds; true labels,
covariances and latent responses do not initialize or select candidate branches.
Two independent R processes use alternating fields with fixed case seeds and
distinct checkpoint paths. Every layout/field/job is separately fingerprinted.

Sixteen attempts have separated candidate likelihoods at the prescribed
numerical resolution; eight remain comparison_unresolved with all overall
parameters NA. The sixteen numerical winners comprise ten positive-range
spatial candidates and six independent-noise endpoints. Thus fourteen overall
range outputs are missing (six endpoints plus eight unresolved). Every one
of the forty-eight endpoint branch range outputs is missing. Branch-specific
means/SDs/nuggets and all candidates remain visible; an unresolved comparison
does not mean that no observations or no location information exist.

Two one-exact cases are private research and are explicitly marked in the plots.
A secondary description using the existing public count rule has twenty-two
eligible attempts, fourteen numerical winners and eight unresolved comparisons.
The primary results retain all cases; no new count filter or eligibility waiver
was introduced. Numerical separation is not scientific validity or evidence
that the selected family generated the observations.

Most importantly, seven endpoint-generated analyses still return a positive
range numerical winner: six under site-only generation and one under independent
noise. The other endpoint-generator outcomes are five independent-noise winners
and four unresolved. These counts on shared layouts are descriptive, not a
calibrated false-positive rate. Maximizing a more flexible covariance likelihood
can fit sampling variation; a fitted positive range cannot be relabeled as
detected cross-site/geological signal. No site-only candidate is the unique
numerical winner in this small study, despite correct checked block likelihoods;
finite kernels can approximate that endpoint and fit incidental cross-site
patterns. A separate statistical evidence procedure is needed.

There are 36,346 optimization evaluations: 33,935 finite and 2,411 explicit
outside-search trials. No within-search likelihood failure occurred in this
pilot and no finite probability floor is used. Eight new 32,768-draw batches
score frozen branch parameters. Site/independent endpoint scores are deterministic;
maximum observed independent-engine discrepancy is 1.06e-11. Numerical comparison
margins combine paired integration MCSE, deterministic disagreement and a small
numerical floor. They are resolution heuristics, not hypothesis tests, covariance
selection guarantees or confidence cutoffs. Deterministic integration MCSE zero
does not imply zero sampling uncertainty.

Evidence:

- [Endpoint findings and statistical limits](output/spatial_endpoints_fresh_2026-11-01/INTERPRETATION.md)
- [Full report and figures](output/spatial_endpoints_fresh_2026-11-01/README.md)
- [Overall outputs](output/spatial_endpoints_fresh_2026-11-01/endpoint_results.csv)
- [All seventy-two branch candidates](output/spatial_endpoints_fresh_2026-11-01/branch_candidates.csv)
- [Every range value or NA and its reason](output/spatial_endpoints_fresh_2026-11-01/range_reporting_audit.csv)
- [Numerical comparison gaps/margins](output/spatial_endpoints_fresh_2026-11-01/numerical_comparisons.csv)
- [Pre-validation numerical development](output/spatial_endpoints_fresh_2026-11-01/NUMERICAL_DEVELOPMENT.md)

Overall CSV reconstructs byte-for-byte from all twenty-four saved case rows.
All twelve generating covariance matrices, paired censoring and two-record
coordinate supports were independently checked. Completed init/worker1/worker2/
finish entry points preserve all 76 then-existing artifacts byte-for-byte.
Source/protocol/runtime and design hashes match; source snapshots and final
SHA-256 artifact manifests are retained. Both final figures were visually
inspected, with private one-exact cases marked. The initial report source/README
remain preserved. Only excluded research scripts and development records changed;
no public package function, confidence limit, count rule, dependency/version,
real USGS estimate, release, commit or push changed. Broader package checks were
not redundantly repeated for excluded research code, and old cloud-only archives
were not reread. Original data, papers and all prior numerical failures remain
intact. Confidence limits remain experimental.

Next gate: design and independently validate statistical evidence for cross-site
dependence under site-only/independent null settings, preserving site blocks,
censoring and survey geometry. The numerical gap criterion is not that test.
Any procedure informed by this pilot requires another frozen fresh study; these
responses cannot be its independent confirmation set. No automated signal claim
or population confidence interval is introduced from these point-estimation
and numerical checks.


## Site-block spatial-order evidence initiated — November 1, 2026

Following the endpoint pilot, implemented a separate research test for positive
spatial ordering under an explicitly exchangeable whole-site null. It uses only
comparisons determined by observed intervals, never substituted concentrations
or fitted covariance parameters. Integer dominance scores sum definite between-
site record orderings within each paired site. Permuting entire site blocks
preserves repeats, observed bounds and common censoring design; geometry stays
fixed. The statistic is the maximum Moran-type quadratic over three fixed
scales, and that same maximum is recomputed on every permutation. The p-value
includes the observed identity with conservative ties and plus-one counting;
it cannot be zero. This is not a spatial bootstrap confidence procedure.

Reporting metadata for every record and an explicit exchangeability-under-null
assumption are required. Missing bounds, heterogeneous reporting designs and
absent score variation return NA rather than deleting observations or filling
values. Correlated repeats are allowed under the site-only null; heterogeneous
geology/methods, trends or response-selected locations can invalidate inference.
No real USGS p-value is computed in this study. Informative censoring bounds
can support this test without a continuous concentration MLE; existing public
concentration-fit count rules are unchanged.

Thirteen rank/equivariance/exhaustive-level/plus-one numerical comparisons and
missing/design/RNG guards passed before validation responses. The exact six-site
permutation group (720 assignments) was checked for conservative rejection at
four levels, including ties and multiple distance scales. Mathematical primary
sources consulted: Hemerik & Goeman (2018), doi:10.1007/s11749-017-0571-1, and
Phipson & Smyth (2010), doi:10.2202/1544-6115.1585. These justify the randomization
machinery under invariance, not our new score's power or real-NGDB eligibility.
This differs from the earlier Kelly warning: here permutation deliberately
calibrates an exchangeable null, rather than pretending to preserve spatial
dependence for a confidence/association placebo exercise.

Frozen protocol: output/spatial_order_signal_2026-11-01/PROTOCOL.md. Pre-generated
400 fresh layouts/plans and 2,000 fields before any outcomes. Development has
50 layouts/regime and confirmation 150/regime, with uniform/clustered geometry,
24 sites/two records per site. Scenarios: independent, site-only, two positive
spatial ranges and a deliberately heterogeneous-reporting refusal check; paired
moderate/heavy censoring produces 4,000 attempts. No tuning after development is
allowed; confirmation requires identical frozen source. Test unavailable values
stay NA, with availability separate from operational rejection summaries. No
production API, covariance-family declaration, data fill or CI is introduced.


## Site-block spatial-order evidence completed — November 1, 2026

Completed the frozen development and independent confirmation phases without
method changes or tuning. The study contains 400 fresh layouts, 2,000 response
fields and 4,000 paired-severity attempts. Development used 50 layouts per
geometry regime and confirmation 150; partners within layouts are dependent.
No covariance model was fitted by this evidence test.

Confirmation false-rejection rates across the independent and site-only null
strata were 3.3–7.3% at nominal 5%. Every stratum's 95% Monte Carlo interval
includes 5%: compatible with nominal control at this simulation precision,
not proof of exact empirical calibration or general NGDB applicability. Pooling
four null outcomes within each layout first gives 6.0% for clustered layouts
(layout Monte Carlo SE 1.15 percentage points) and 6.17% for uniform layouts
(SE 1.06 points). Development uniform pooling was 9.0% (SE 2.34 points);
this result is retained and did not trigger retuning or additional selective
sampling. The 600 null tests per confirmation regime are not counted as
600 independent surveys.

Confirmation operational detection rates were 38.7–78.0% with moderate
censoring and 21.3–35.3% with heavy censoring. Non-rejection does not establish
independence or adequate spatial coverage. Across both phases, 3,129 of 3,200
eligible-design attempts were available; 71 had no order-score variation and
kept NA p-values, including 51 confirmation cases. All 800 deliberately
heterogeneous reporting-design cases were refused with NA. Twelve available
tests had zero exact records but informative censoring-bound orderings; this
does not alter public concentration-fit count rules or establish an MLE.
No unsampled concentration, missing reporting limit or missing result was filled.

Validation and retained evidence:

- Thirteen numerical/exhaustive comparisons and missing/design/no-score/RNG
  guards passed before fresh responses, including all 720 six-site assignments.
- All 4,000 saved tests were checked against raw bounds, dominance scores,
  permutation scale maxima and plus-one counts. Independent/site-only null
  covariance zeros were verified for every generated layout.
- Development/confirmation CSVs reconstruct byte-for-byte from all 400
  checkpoints. Completed init, both worker modes and finish preserve all 829
  then-existing study artifacts byte-for-byte.
- Seven frozen source/input and 401 design fingerprints match. Four canonical
  source snapshots, runtime information and final artifact hashes are retained.
- Both figures were visually checked. Percentage-axis labels were clarified
  in reporting only; the frozen test, protocol, fields and results were unchanged.

[Scientific interpretation](output/spatial_order_signal_2026-11-01/INTERPRETATION.md),
[full report and figures](output/spatial_order_signal_2026-11-01/README.md),
[all results/statuses](output/spatial_order_signal_2026-11-01/all_results.csv),
[availability and simulation-rate intervals](output/spatial_order_signal_2026-11-01/rejection_and_availability.csv),
[layout-grouped null summaries](output/spatial_order_signal_2026-11-01/pooled_null_layout_summary.csv).

Primary randomization references consulted online were Hemerik & Goeman (2018)
and Phipson & Smyth (2010); scope and links are recorded in PAPER_REVIEW.md and
PAPER_FETCH_LIST.md. These support the invariant-null randomization machinery,
not real-data eligibility or power of the new score. No new blocking supplement
request was added, and online consultation was not described as receipt of a
local Papers file.

Only excluded research scripts, their outputs and development/research notes
changed. No public R function, dependency, version, release or real-USGS estimate
changed. The previous successful package check remains the latest package
validation; broad checks and reads of old cloud-only archives were not repeated.
Original data, papers and earlier failures remain intact. No commit or push.
Confidence limits remain experimental. This test does not select a covariance
family/range, correct preferential sampling or validate population confidence
limits; across-element/domain multiplicity is outside this single-domain study.

Next gate: a real-data eligibility audit of comparable medium/method/support
domains, complete reporting-limit design metadata even for exact records, and
defensible site grouping/exchangeability. Matching metadata alone does not prove
exchangeability. Where required information is missing, evidence stays NA; do
not infer a reporting limit from the smallest exact value or automatically
apply the test to every USGS element.


## USGS spatial-order eligibility audit — November 1, 2026

Completed a documentation, access and readable-metadata audit rather than
applying the experimental p-value to the mixed USGS compilation. Both bundled
data dictionaries and schema files are readable; all 102 documented columns
across six table schemas agree. The soil metadata header and complete table
were also read. Initial short reads timed out; a later bounded read retrieved
soil main.csv. The two bestvalue tables, two all-chemistry tables and sediment
metadata still timed out in the audit (after earlier bounded attempts), and
remain explicitly unavailable for new record-level analysis. Their contents,
chemistry counts and eligibility are NA, not zero or assumed missing columns.
No broad old cloud-only archive reads were performed.

Readable soil metadata contains 104,181 total sample records and 102,035 U.S.
records, all with valid coordinate pairs. These span 59,368 exact coordinate
groups; 12,559 repeated groups contain 55,226 sample records. Maximum records
at a coordinate pair is 854; 2,526 groups exceed four records, accounting for
32,220 records. These are coordinate matches, not established physical sites,
survey repeats or independent locations. No records were dropped, thinned,
averaged, moved or assigned invented coordinates/site labels.

Repeated-coordinate groups have multiple recorded depth labels in 2,698 cases,
size-fraction labels in 2,491, and sample-zone labels in 1,912. Counts overlap;
raw label differences can be spelling/units as well as substantive differences.
At least one depth is blank in 9,469 repeated groups and a size fraction blank
in 7,068. Missing labels are not agreement. Nonblank soil coverage is 12.64%
for depth, 24.93% for size fraction, 3.85% for preparation and only 21 records
(0.0206%) for collection protocol; 30.97% have a sample-zone entry. Nonblank
values are not certified usable/comparable metadata. The dictionary distinguishes
`sample_zone` (horizon zone) from `horizon` (how well defined the horizon was).

Documentation declares method/unit/species and observed censoring encodings,
but no dedicated lower/upper reporting-design fields for every record or
record-level limit provenance. This does not establish that the information
cannot be recovered from external method/laboratory/era documentation later.
A censored observed bound is distinct from a fully specified design for exact
records. `lab_id` identifies a sample and `job_id` a batch, not a laboratory;
bestvalue `rec_id` is table-local, so no cross-table join is inferred. A
physical-site/repeat map, defensible exchangeability, real metric distance
scales and across-domain multiplicity also remain unestablished. Simulation
unit-square scales are not automatically degrees/metres.

The independent R check matches all soil row totals, reported field counts
and the entire coordinate-multiplicity distribution from the Python audit;
U.S. metadata sample IDs are unique. Its initial verification attempt stopped
before output on a data.table column/loop-variable name collision. The failed
source and correction note are preserved; corrected variable names passed
reconciliation without changing primary counts. Both current sources and all
readable documentation/metadata fingerprints are retained. No unavailable
input was hashed as if it had been read.

[Audit findings and scope](output/usgs_spatial_order_eligibility_2026-11-01/README.md),
[thirty eligibility requirements](output/usgs_spatial_order_eligibility_2026-11-01/eligibility_requirements.csv),
[input access status](output/usgs_spatial_order_eligibility_2026-11-01/input_availability.csv),
[repeated-coordinate support review](output/usgs_spatial_order_eligibility_2026-11-01/repeated_coordinate_support_review.csv),
[independent reconciliation](output/usgs_spatial_order_eligibility_2026-11-01/independent_metadata_verification.csv),
[next data requirements](output/usgs_spatial_order_eligibility_2026-11-01/DATA_REQUIREMENTS.md).

All real-data p-values and the number of eligible domains stay NA. This is a
completed access/documentation/soil-metadata audit, not a completed all-element
chemistry-domain review, a zero-eligibility conclusion, or evidence of absent
spatial structure. No new numerical threshold, preferred domain or concentration
estimate was introduced. Package functions/dependencies/version and previous
USGS results remain unchanged; only excluded research and development records
changed. Existing package-check results remain current; no release, commit or
push. Original data and papers are read only. Confidence limits remain
experimental; no additional paper/supplement is required for this audit.

Next work: when record files are readable, inventory method/unit/project/support
domains without selecting on the proposed test statistic. Establish traceable
reporting designs and physical-site definitions before any real-data test.
Where the needed information remains unavailable, preserve NA and the reason.


## V1.1 release candidate prepared — November 1, 2026

Froze release scope around USGS censoring/data integrity, multi-element workflows,
cell-declustering diagnostics, information/support review and empirical censoring
bounds. Spatial covariance/order/HAC methods remain excluded research; the
opt-in block-bootstrap API is explicitly experimental. Routine weighted limits
remain NA, and ordinary model limits require independent observations and an
appropriate lognormal model. No spatial confidence procedure was promoted.

Release review found the older exported `fit_lognorm_with_quants()` still had
a separate interval construction, a minimum-left-censor rule and insufficient
failure handling. It now uses the same interval validation and shared numerical
fit guard as the main estimator, supports uncensored/right-only samples, honors
custom probabilities, always reports the model median, and returns explicit
insufficient/failed statuses with all unavailable fitted quantities missing.
Its successful default column names/order and mg/kg concentration units remain
unchanged. This safety/contract correction is described in NEWS and documented;
new regression tests check analytic uncensored MLEs, units, custom quantiles,
invalid limits/conflicting flags, count guards and warning suppression. Model
independence assumptions are explicit in both public estimator help pages.

Updated source/citation metadata to 1.1.0 and rewrote NEWS as release notes,
including compatibility changes for duplicate widening and quantile statuses.
CITATION's previous date-released was removed rather than claiming an unpublished
candidate had been released. README and release checklist identify the candidate.
Roxygen documentation was regenerated. The existing five-job CI workflow was
compared with the official [r-lib standard template](https://github.com/r-lib/actions/blob/v2/examples/check-standard.yaml).
No action/dependency upgrade was needed for this preparation.

The clean source stage has 97 text files and the built archive 101 files,
including generated vignette material. The installed archive excludes PFAS,
sibling packages, papers, development research, data and generated outputs.
The repository publication allowlist separately includes the research source
and development notes, with their local/unbundled artifact limits documented.
Unrelated workspace documents and sibling packages remain untouched. The first
build command used an incorrect relative directory; the correct invocation
then built successfully. Review found a local `.DS_Store` in that initial source
inventory (R CMD build omitted it); preserved the initial stage/archive and
rebuilt from a text-only allowlist before the final check. A credential-pattern
scan of allowlisted text found no matches; this is not a general security audit.

Final local tests: 77 tests, 410 passing expectations, zero failures/errors/
warnings/skips. `R CMD build --no-manual` and `R CMD check --no-manual` both pass,
Status: OK, with vignette creation and rebuilding. Installed 1.1.0 single- and
all-element example scripts ran on synthetic chemistry/metadata. Leading-zero
IDs survive; all-left/no-exact Au and unsupported-unit Zn remain present with
unavailable estimates; all six declustering rows retain missing routine limits.
The single-element plot was visually checked. No original NGDB estimate was
refit during this release gate.

[Release checklist](RELEASE_CHECKLIST.md),
[candidate check log](output/release_candidate_1.1.0_2026-11-01/check.log),
[test summary](output/release_candidate_1.1.0_2026-11-01/test_summary.csv),
[publication allowlist](output/release_candidate_1.1.0_2026-11-01/publication_allowlist.csv),
[installed-example verification](output/release_candidate_1.1.0_2026-11-01/CLI_VERIFICATION.txt).

The GitHub connector was unavailable. The restricted CLI initially reported
invalid authentication; the same read-only authentication check through the
approved network connection succeeded, so credentials were not reset. Remote
main matches starting commit 911f8aadb294b825ed492eb5393f1a71688c8462 and the
intended codex/release-1.1.0 branch did not yet exist when checked. Publication
will use the verified allowlist and a review branch, with remote checks before
any tag. No merge/tag/release is authorized by this preparation record.


## V1.1 review branch pushed and matrix passed — November 1, 2026

Committed only the explicit 180-file publication allowlist and pushed
`codex/release-1.1.0`. Source candidate commit is
`c7845c3bc1568f2078a9e135d0da0e756113c8ca`; all allowlisted committed blobs were
verified byte-for-byte against the publication manifest before pushing. The
commit changes 140 files; unchanged allowlisted files are also verified.
Unrelated untracked files and original data remain local and intact. Remote
main was unchanged at publication; no merge, tag or release was created.

Opened [draft PR #1](https://github.com/diskoisland/geoCensor/pull/1).
All five jobs in [the source-candidate workflow](https://github.com/diskoisland/geoCensor/actions/runs/41385826608)
passed: macOS/current R, Windows/current R, Ubuntu/current R, Ubuntu/oldrel-1
and Ubuntu/devel. The devel job took 24m31s, chiefly dependency setup (about
23 minutes), then passed its package check. This delay was not a failed check.
Saved the complete job/step result in the local release-validation artifacts.
The initial commit appeared delayed in the cloud-backed checkout, completed
successfully, and was verified before push; no failed commit or data reset was
used. No GitHub credential change was needed.

A final documentation-only commit records the actual publication and matrix
outcome in this log and RELEASE_CHECKLIST.md. All 97 checked package source
files and the checked archive SHA-256 remain identical. Head checks must still
be green before tagging; this record is not a release announcement. The draft
is ready for review, and spatial confidence limits remain experimental.


## V1.1 publication finalization — November 1, 2026

The user authorized moving from the reviewed candidate to merge and V1.1 tag/
release publication. Fresh inspection found PR #1 open, draft, mergeable and
clean, with all five checks passing on head 1a78cc671a9f70ededb65d7533ce1db0778a0bf0.
The GitHub connector remained unavailable; the signed-in CLI is the fallback.

Removed candidate/unpublished wording from README and NEWS, added the client-
dated November 1, 2026 citation release date, and refreshed this log/checklist.
No public function, dependency, test, example, vignette or scientific inference
rule changed. Ordinary model assumptions and experimental spatial limits remain
explicit. The final installed-source changes are only README/NEWS wording.
Fresh final archive build/check passes, Status: OK, including package tests
and vignette creation/rebuilding. The final stage has 97 files and archive 101;
PFAS, data, papers, private research scripts and outputs are excluded.

The publication procedure rechecks the exact final PR head, preserves branch
protections, performs a normal merge without force/deletion, checks merged-tree
and file hashes, and waits for all five main-branch checks before creating an
annotated v1.1.0 tag. It publishes only the checked archive and traceable
validation/checksum assets, then verifies downloaded archive bytes. Final head,
merge commit, workflows, tag and asset verification are retained in the release's
VALIDATION.md and local publication receipt. Existing research failures and
original data remain intact; spatial confidence limits remain experimental.

Final archive SHA-256: `bfd47d1d89828b8cd87efe04a835cdbddb241cefaa12112d9e919fb1cb1232e2`.

[Release publication record](https://github.com/diskoisland/geoCensor/releases/tag/v1.1.0)
[Final local check](output/release_1.1.0_2026-11-01/check.log)
[Publication receipt](output/release_1.1.0_2026-11-01/PUBLICATION_RECEIPT.md)
