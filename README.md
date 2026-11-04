# geoCensor

**Censoring-aware analysis of USGS National Geochemical Database (NGDB) data.**

Development source version **1.2.0.9006**. Published versions and downloadable assets are listed
on the [releases page](https://github.com/diskoisland/geoCensor/releases). See
[NEWS.md](NEWS.md) for changes and compatibility notes.

I wrote geoCensor for working with USGS National Geochemical Database tables where a lot of the useful information is tucked inside detection limits, qualifiers, and non-detects. Those details matter: in the NGDB, detection limits vary by method and over time, so dropping non-detects or replacing everything with LOD/2 can move the center of the data in ways that are hard to see later.

The package turns USGS `qvalue` and `qualifier` fields into concentration intervals, including open-ended intervals for non-detects and greater-than values. It then estimates geometric means and confidence intervals with maximum-likelihood lognormal models using interval-censored survival regression. It also includes LOD/2-style substitution summaries for comparison, distribution diagnostics, spatial joins, robust correlation analysis, robust PCA, and mapping helpers.

To generate the MLE:
1. non-detects were coded as `(0, LOD]`;
2. greater-than values were coded as `(ULOD, Inf)`; and
3. unqualified detections were treated as point observations under the working model.

That last coding rule does not certify measurement precision. Documented
semiquantitative methods retain review flags; reporting-class bounds are not
inferred from a method label.

The lognormal model is fit with `survival::survreg(..., dist = "lognormal")`, and `exp(meanlog)` is used as the geometric mean.

Ordinary model confidence limits assume independent observations and an
appropriate lognormal model. They do not account for spatial dependence or
preferential sampling. Declustering-weighted fits omit routine confidence
limits; spatial confidence procedures remain experimental.

The censored-data treatment follows Helsel (1990), which contrasts fixed
substitution with distributional maximum-likelihood and robust methods:
D. R. Helsel, “Less than obvious: Statistical treatment of data below the
detection limit,” *Environmental Science & Technology* 24(12), 1766-1774,
[doi:10.1021/es00082a001](https://doi.org/10.1021/es00082a001).


A runnable, synthetic end-to-end example is available in the
[USGS workflow vignette](vignettes/usgs-workflow.Rmd). Installed packages can
open it with `vignette("usgs-workflow", package = "geoCensor")`.

## Installation

Release changes and compatibility notes are recorded in [NEWS.md](NEWS.md).
The command below installs the current public main branch. Check the installed
version before using the versioned example commands. A release candidate can
also be installed directly from its checked source archive.

```r
# install.packages("remotes")
remotes::install_github("diskoisland/geoCensor")
packageVersion("geoCensor")

# Alternatively, after obtaining the checked 1.2.0 source archive:
# install.packages("geoCensor_1.2.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(geoCensor)

# 1. Parse a chemistry table with lab_id, species, unit, qvalue,
#    and qualifier columns into censoring fields for one element,
#    converting supported units to ug/kg
ag <- build_ag_from_csv("ngdbsoil/bestvalue.csv", element = "Ag", to_ugkg = TRUE)

# 2. Attach sample metadata (state, coordinates) by lab_id
main <- read_clean_main("ngdbsoil/main.csv")
ag_main <- join_soil_main(ag, main)

# 3. Censoring-aware lognormal MLE geometric mean, by state
gm_by_state <- mle_by_state(ag_main)
head(gm_by_state)        # state, n_total, n_left, n_right, geomean, ci_low, ci_high

# 4. Spatial coverage diagnostics for the same state metrics
coverage_by_state <- spatial_coverage_by_state(ag_main)
gm_by_state <- dplyr::left_join(gm_by_state, coverage_by_state, by = "state")

# 5. Descriptive pooled fit for the supplied records; review the domain first
pooled <- fit_cen_lognorm(ag_main)
pooled

# 6. Inspect the fitted model when estimation is available
if (pooled$method == "survreg_lognormal_MLE_L+R") {
  plot_lognormal_fit(ag_main)
}

# 7. Map / rank
plot_ranked_bar(gm_by_state)
plot_us_map(gm_by_state)
```

## Example scripts

Installed example scripts show how to build Ag MLE state tables and hotspot
maps from NGDB `bestvalue.csv` files:

```r
system.file("examples", "build_ag_mle_state_df.R", package = "geoCensor")
system.file("examples", "ag_seed_state_mle.R", package = "geoCensor")
system.file("examples", "test_ag_spatial_coverage_step_by_step.R", package = "geoCensor")
system.file("examples", "ag_state_mle_spatial_coverage.R", package = "geoCensor")
system.file("examples", "test_ag_lognormal_fit_plot_step_by_step.R", package = "geoCensor")
system.file("examples", "test_idaho_ag_lognormal_fit_plot.R", package = "geoCensor")
system.file("examples", "test_idaho_ag_weibull_fit_plot.R", package = "geoCensor")
system.file("examples", "test_us_pooled_ag_distribution_fit_plot.R", package = "geoCensor")
system.file("examples", "test_us_pooled_ag_nonparametric_cdf_plot.R", package = "geoCensor")
system.file("examples", "ag_lognormal_fit_plot.R", package = "geoCensor")
system.file("examples", "ag_hotspot_map.R", package = "geoCensor")
```

For a single element, run the reusable command-line workflow against an NGDB
chemistry or best-value CSV:

```sh
Rscript inst/examples/analyze_element_dataset.R data/ngdbsoil/bestvalue.csv Ag output/Ag
```

This requires installed geoCensor 1.2.0.9006 or newer. Replace `Ag` with another element
symbol. The workflow converts supported units (`ppm`, `mg/kg`, `ppb`, `ug/kg`, micro-symbol equivalents, and percent units) to ug/kg and
exports normalized observations, the pooled estimate, censored distribution
comparisons, a fitted CDF plot when estimation succeeds, serialized results,
a validation summary with exclusion reasons, a run summary, and R session information. Unsupported units stop the analysis.

Input validation retains invalid chemistry rows with missing concentration fields
and an explicit `validation_status`: missing or invalid values, nonpositive
detections, invalid limits, conflicting encodings, unknown qualifiers, or
unsupported units. Blank/missing, `L`, `N`, `<`, `G`, and `>` qualifiers are accepted;
other codes need review before use. Read chemistry with `read_usgs_chemistry()`
or `build_ag_from_csv()`. Both installed analysis commands use the same reader.
All chemistry columns remain character, preserving leading-zero/large identifiers,
numeric text and sparse qualifiers. Blank fields and literal NA are missing;
other text is retained. Malformed CSV records and duplicate/empty headers stop
before analysis. A value such as `4.11110` retains its trailing zero and cannot
acquire an upper-censoring suffix through numeric type guessing.

```r
chemistry <- read_usgs_chemistry("ngdbsoil/bestvalue.csv")
```

The all-element command reads metadata as text and parses only explicitly
supplied projected-coordinate columns as numeric. Blank/NA and explicit nonfinite
coordinates remain unlocated. Invalid coordinate text stops
before outputs are saved. Supply an appropriate CRS and coordinate units; neither
is inferred. Source-precision exports retain character source values/IDs rather
than guessed numeric columns. In-memory builder/analysis APIs retain their
existing input-type behavior.

Metadata joins reject duplicate or blank metadata IDs, warn about unmatched
or repeated chemistry IDs, and retain chemistry rows without multiplying them.
The `join_audit` attribute records match counts. Wide metals tables reject
repeated sample/element pairs instead of averaging them; choose the appropriate
result explicitly before widening. Their `validation_summary` attribute records
usable and excluded input counts.

## Save a reproducible metadata audit

```r
run <- audit_sampling_file(
  "metadata.csv", "output/support_review", group_cols = "state",
  x_col = "long_wgs84", y_col = "lat_wgs84", include_memberships = TRUE
)
```

This workflow requires1.2.0.9005 and a new or empty output directory. It exports
support/membership tables, raw metadata, coordinate conversions, RDS, input MD5,
configuration and session information. Fields are read as text to retain IDs and
formatting. The analysis copy treats trimmed blank/NA as missing; exported raw
metadata and membership source columns preserve original text. Location/group
coordinates are numeric, with each conversion recorded separately.

The installed `examples/audit_sampling_support.R` command provides the same run:

```text
Rscript audit_sampling_support.R metadata.csv output/support_review --group-cols=state --coords=long_wgs84,lat_wgs84 --memberships
```

Use `system.file("examples", "audit_sampling_support.R", package = "geoCensor")`
to locate the installed script. Optional `--support-cols=depth,mesh_size` and
`--id-col=lab_id` select fields explicitly. Without coordinates or memberships,
the command exports the ordinary field/combination audit. Invalid coordinate
text stops before saving; blank/NA and nonfinite coordinates remain unlocated.
Existing nonempty destinations and duplicate/unknown options are rejected.
No CRS, units, dates, physical identity, independence or dilution is inferred.

## Link records to recorded-position review groups

```r
support <- sampling_support_audit(
  metadata, group_cols = "state", x_col = "long_wgs84", y_col = "lat_wgs84",
  include_memberships = TRUE
)
support$memberships
support$coordinate_groups
```

This opt-in requires1.2.0.9004. It preserves every original metadata row/column,
including raw IDs and unlocated records. Exact numeric coordinate pairs define
global review groups across analysis domains. Existing domain-specific locations
and the other default support-audit components remain unchanged. Coordinates must
use a common reference system; no transformation or physical-site identity is inferred.

Missing coordinates remain unlocated with NA group/count fields. Position labels
are cohort-specific, stable under row reordering, and are not permanent site IDs.
Repeated positions receive a review flag; singletons are not certified independent.
No records are collapsed or assigned weights. These memberships can support
explicit decisions about which records to hold together in later validation.

## Audit source-summary precision

For best-value tables with a USGS `summary` column, compare the reported number
with its first explicitly associated source entry:

```r
source_review <- audit_summary_precision(chemistry)
source_review[c("lab_id", "species", "unit", "qvalue",
                "summary_precision_status", "summary_numeric_delta")]

reviewed_run <- analyze_elements(chemistry, audit_source_precision = TRUE)
reviewed_run$source_precision_summary
```

Batch auditing requires geoCensor **1.2.0.9002** or newer; the current installed
commands require **1.2.0.9006** so they use guarded recovery and chemistry imports.
The opt-in adds `source_precision` and `source_precision_summary`; existing
estimates, validation, method tables and spatial weights stay unchanged. The
row table retains requested input order, raw IDs/values/units and all original
chemistry columns. Metadata-only grouping columns are appended when requested.
Candidates and numeric differences retain each row's original units, while
estimates continue to use ug/kg. Counts include input rows excluded from fitting.

Differences are source-representation evidence, not estimated analytical error.
No candidate replaces an estimator input. Missing summaries, malformed entries,
method mismatches, UNKNOWN labels and coding conflicts stay explicit. Existing
enabling the audit does not change parser recovery; it compares the original
input before recovery. A matching value does not certify analytical accuracy.

Both installed commands accept `--audit-source-precision`, alone or with
`--stratify-methods`; the all-element command can also combine it with
`--method-overlap`. They export `source_precision.csv` and
`source_precision_summary.csv` plus the corresponding RDS components only when
requested. Choose a new output directory for each run.

## Review analytical methods

```r
reviewed <- flag_analytical_methods(chemistry)
reviewed[c("amethod", "method_review_status", "review_semiquantitative")]
```

The USGS codes `ES_SQ` and `ES_H2O_P` are documented semiquantitative methods.
They receive a review flag and documentation source. Missing/blank methods,
`UNKNOWN`, and other codes remain unresolved; an NA flag does not certify a
quantitative measurement. Original values, method labels, censoring fields and
validation statuses are preserved. The helper supplies no numeric reporting-class
bounds, detection limits or laboratory/era assumptions.

`build_element()` adds these fields automatically. `analyze_elements()` adds
method counts/status to estimates and a `methods` table with input and unweighted
fit-subset counts. Declustering summaries use their own chemistry/coordinate
subset, so a flagged record without usable coordinates is not counted as used
in a spatial fit. Installed single/all-element commands export `methods.csv`.
Direct `fit_cen_lognorm()` and `fit_lognorm_with_quants()` calls, and
`mle_by_state()`, also retain method counts/status. Custom censoring columns
and zero-weight exclusions use the same eligible subset as their fit.
“Used” means eligible input even when the fit is insufficient or fails.
Plots include method-review counts in saved captions and reject numerical
fit failures before drawing a fitted curve. When every usable interval includes
the same uncensored-coded value, a positive-scale lognormal maximum does not
exist: the fitted mean, limits and quantiles remain missing, and no fitted curve
is drawn. Records, method flags and empirical censoring bounds remain available.
Distribution comparisons retain method counts/status for every candidate,
including failed or insufficient fits. Empirical CDF bounds append the same
review fields to `support` using their own valid, positive-weight subset.
CDF plots retain these cautions; older saved results explicitly report missing
method-review information rather than reconstructing it.

Count-based review remains separate: `n_exact` counts uncensored input coding,
not verified measurement precision. A result can clear the count thresholds
and still require method review. No estimates are filtered or changed by these
annotations, and spatial confidence limits remain experimental.

### Compare recorded methods before pooling

```r
method_check <- analyze_elements(
  chemistry, main_df = metadata, group_cols = "state",
  stratify_methods = TRUE
)
method_check$method_estimates[c("element", "state", "method_code", "n_input",
  "n_total", "geomean", "method_information", "review_status")]
```

This optional table separates unweighted summaries by normalized method code
within each element/domain. Matching trims spaces and ignores case; original
labels stay intact. Blank/absent methods remain an unresolved NA bucket, with
`UNKNOWN` separate. A missing-method bucket does not establish a shared
procedure. Unusable inputs, insufficient groups and failed fits remain explicit.
The six existing result tables and pooled estimates are unchanged. Method
summaries use chemistry eligibility, including records without usable coordinates.

Differences between method strata can also reflect location, geology, sample
support, era or preferential sampling. They are descriptive sensitivity results,
not a calibration or an estimate of analytical bias. Ordinary confidence limits
still assume independent observations and a suitable lognormal model.

Both installed commands accept `--stratify-methods` and then export
`method_estimates.csv` as well as the usual outputs. This feature requires
version 1.2.0.9006 or newer.
Choose a new or empty output directory for each command run; existing results
are preserved so optional tables and plots cannot be confused with a later run.

### Check whether method groups cover the same area

For one selected element, medium and domain, use valid chemistry with projected
coordinates to inspect overlap at several cell widths:

```r
overlap <- method_spatial_overlap(
  analysis_data, cell_sizes = c(5000, 10000, 20000),
  x_col = "x_m", y_col = "y_m", n_offsets = 2
)
overlap$support
overlap$methods
overlap$pairs
```

All methods share the same grid at each placement. Pair summaries report shared
cells, their union and the shared/union fraction (Jaccard), with minima/means/
maxima across placements. These ranges describe grid sensitivity, not confidence
limits. Missing methods remain explicit buckets outside named-method pairs.
When a code lacks usable spatial data, comparisons stay NA rather than becoming
zero overlap. Genuine zero overlap is reported only when both codes have data.

Cells shared by two method labels do not establish paired samples, comparable
physical support or analytical agreement. Exact coordinate pairs may be rounded
or reused and are not independent-site counts. This descriptive audit fits no
concentration model and does not resolve preferential sampling or unsampled areas.

To run the same audit for every selected element and domain:

```r
geographic <- analyze_elements(
  chemistry, main_df = metadata, group_cols = "state",
  cell_sizes = c(5000, 10000, 20000), x_col = "x_m", y_col = "y_m",
  n_offsets = 2, audit_method_overlap = TRUE
)
geographic$method_overlap_support
geographic$method_coverage
geographic$method_overlap
```

This appends three labelled audit tables and leaves existing outputs unchanged.
It is independent of `stratify_methods = TRUE`; both options can be combined.
The audit uses the supplied `n_offsets` (batch default 10, or 100 placements),
whereas the standalone audit defaults to 2. Review cell sizes and placement counts
for the intended spatial domain. The batch audit uses sorted unique cell sizes.

The all-element command accepts `--method-overlap` when metadata, projected
coordinate columns and cell sizes are supplied. It exports
`method_overlap_support.csv`, `method_coverage.csv` and `method_overlap.csv`, plus
the corresponding RDS components. This command now requires version 1.2.0.9006
or newer. Missing method/coordinate information stays explicit in those tables.

## Audit sampling support before pooling

```r
support <- sampling_support_audit(
  metadata, group_cols = "state",
  support_cols = c("depth", "sample_zone", "horizon", "mesh_size"),
  x_col = "long_wgs84", y_col = "lat_wgs84"
)
support$fields
support$locations
```

Use fields actually available in your dataset; `support_cols = NULL` selects
available standard metadata fields. In USGS soil exports `sample_zone` records
the horizon zone; `horizon` records how well defined it is. Sediment metadata
may omit both. This
is a descriptive audit: different recorded labels and missing information are
reported separately. Exact coordinate reuse may reflect rounding or different
sample supports. Matching or missing labels do not establish comparability;
no records are automatically merged, dropped, or assigned weights.

## Analyze all elements

Analyze every nonblank species in an NGDB chemistry or best-value table:

```r
results <- analyze_elements(chemistry)
results$estimates
results$validation
results$methods
results$distributions
```

Batch estimates now include `n_exact`, `exact_fraction`, and configurable
count-based review flags. Defaults flag fewer than 10 exact observations or
less than 5% exact observations. These are review conventions, not validated
requirements for reliable inference; clearing the flags does not establish
spatial independence, comparable sample support, or an adequate model.
Use `min_exact_review` and `min_exact_fraction_review` to state alternative
review conventions. Existing estimate tables can be annotated without changing
their estimates with `flag_estimate_information()`.

Use `elements = c("Ag", "Pb", "Zn")` to select a subset. Every selected element
is analyzed separately in ug/kg; unusable units and invalid observations are
counted by reason. Elements with insufficient usable data remain in the results.
Duplicate sample IDs within an element must be resolved before analysis.

For spatial sensitivity, supply metadata with projected coordinates and specify
analysis domains explicitly. Weights are recalculated within each element and
domain, using its own usable observations:

```r
results <- analyze_elements(
  chemistry, main_df = metadata_projected, group_cols = "state",
  cell_sizes = c(5000, 10000, 20000), x_col = "x_m", y_col = "y_m"
)
results$declustering
results$joins
```

Run soil and sediment separately, or supply an appropriate medium column among
`group_cols`. A state is an administrative grouping; use geological domains
instead where they better describe the population being estimated.

The installed package also provides a batch command-line workflow:

```sh
Rscript inst/examples/analyze_all_elements.R bestvalue.csv output/all_elements
# With projected metadata, cell widths in metres, and separate state analyses:
Rscript inst/examples/analyze_all_elements.R bestvalue.csv output/all_spatial main_projected.csv x_m y_m 5000,10000,20000 state
```

It exports estimates, distribution comparisons, validation counts, optional spatial
sensitivity and join audits, serialized results, a run summary, and session details.

## Reducing the influence of clustered samples

`spatial_coverage_by_state()` describes sampling spread. To assess how clustering
affects estimated concentrations, use cell declustering within one sampling
medium and a defined spatial or geological domain. Project coordinates into a
suitable coordinate system first; decimal degrees do not define equal-area cells.
For example, with projected coordinates `x_m` and `y_m` in metres:

```r
sensitivity <- declustering_sensitivity(
  ag_main, cell_sizes = c(5000, 10000, 20000),
  x_col = "x_m", y_col = "y_m"
)
sensitivity
```

Each occupied cell receives equal total influence, shared among its samples.
Weights are averaged over 100 deterministic grid placements. The table compares
weighted censored lognormal estimates against an unweighted estimate on the same
usable chemistry/coordinate subset. Invalid chemistry and missing coordinates
are excluded from both sides of that comparison. Choose cell sizes using sampling
spacing and geological context; inspect several scales rather than selecting a
size solely because it gives the lowest concentration.

The table also compares the raw `exact_fraction` with `weighted_exact_fraction`.
`weighted_left_fraction` and `weighted_right_fraction` complete the shares of
declustering weight across observation types. Many exact records in a dense
cluster may carry little total weight. When almost all weight is censored,
the fitted mean and spread can be sensitive despite successful convergence.
These fractions describe the usable subset; they are not reliability tests or
counts of independent information. They are missing when no rows are usable.

For direct control, calculate weights on your usable analysis subset and pass
them to the estimator:

```r
analysis_data$spatial_weight <- cell_declustering_weights(
  analysis_data, x_col = "x_m", y_col = "y_m", cell_size = 10000
)
fit_cen_lognorm(analysis_data, weight_col = "spatial_weight")
```

The weighted estimator retains left- and right-censoring intervals. Its confidence
limits are intentionally missing: spatial uncertainty needs a separate procedure,
such as a validated spatial block bootstrap. `weight_effective_n` describes weight
imbalance and is not a count of independent locations. Declustering reduces the
influence of densely sampled areas; it does not establish concentrations in
unsampled terrain or remove every source of preferential sampling bias.

## Inspecting distribution shape with censoring

`censored_cdf_bounds()` describes how much of the empirical cumulative
distribution remains unresolved by detection limits. It preserves both censoring
directions and accepts optional observation weights. For one reviewed sampling
medium/domain, with concentration thresholds in the same units as the data:

```r
shape <- censored_cdf_bounds(
  analysis_data, thresholds = 10^seq(-1, 4, length.out = 100),
  weight_col = "spatial_weight"
)
shape$support
plot_censored_cdf_bounds(shape, unit_label = "ug/kg")
```

At each threshold, the plot shows the smallest and largest cumulative fractions
compatible with the observed values and limits. It fits no model and replaces
no censored value. Bounds can be calculated even when a concentration fit is
unavailable. Empty usable subsets have missing bounds.

These describe the supplied observations; they are not confidence bands for a
population. Coincident bounds in uncensored input do not establish measurement
precision or remove sampling and spatial uncertainty. These bounds condition on
the reported censoring codes; semiquantitative class ranges are not inferred.
Inspect distributions within sensible domains and compare
the fitted/nonparametric CDF view from
`plot_lognormal_fit(..., show_nonparametric = TRUE)` where fitting is appropriate.
A model CDF outside an empirical bound is not an automatic rejection rule.
Small-sample non-rejection and a relative AIC winner also do not establish model
adequacy. No ordinary iid normality-test p-value is assigned to these weighted,
censored records.

## Experimental spatial uncertainty

`spatial_block_bootstrap()` is a prototype for assessing uncertainty within one
medium and spatial/geological domain. It samples occupied square blocks with
replacement, recalculates declustering weights, preserves each block's draw
multiplicity, and refits the censored model. Declustering cell width and bootstrap
block width are separate choices. The following widths are illustrative:

```r
boot <- spatial_block_bootstrap(
  analysis_data, x_col = "x_m", y_col = "y_m",
  cell_size = 10000, block_size = 50000,
  n_boot = 999, seed = 42
)
boot$summary
boot$replicates
```

The returned percentile limits have status `experimental_interval`. They are
not automatically added to existing estimates or all-element workflows.
Too few occupied blocks, a failed point fit, or excessive failed bootstrap fits
suppress the limits. Failed replicates remain visible. A seeded call restores
the caller's random-number state.

Block size and origin need sensitivity checks; enough blocks and successful
fits do not establish correct interval coverage. The procedure conditions on
sampled blocks and cannot correct all preferential sampling mechanisms or infer
unsampled terrain. Irregular sampling and correlation across block boundaries
can invalidate simple spatial bootstrap approaches; see
[Lahiri and Zhu (2006)](https://doi.org/10.1214/009053606000000551).
The pilot study in `scripts/analysis/validate_spatial_bootstrap.R` completed 100
datasets in each of five settings with 199 resamples per dataset. Observed
coverage of nominal 95% intervals was 88%, 91%, and 94% for the three
independent-site density settings, and 73% versus 82% for correlated sites with
smaller versus larger blocks. All fits converged. These results show that the
prototype is not calibrated for routine use; it remains experimental and is not
applied to the full USGS results.

A subsequent benchmark study reproduced exact site-log-mean resampling and used
independent confirmation seeds. With 3,000 uncensored confirmation datasets,
analytic t and studentized intervals achieved about 95% coverage; the current
percentile method achieved 94.6%. Censored independent-site coverage was 96–96.5%
in 200 datasets per setting, while correlated-site coverage remained 77–85%.
Basic and normal log-scale formulas did not resolve the spatial shortfall.
The main unresolved issue in those scenarios was underestimated uncertainty from
dependence across blocks. These findings do not validate routine intervals for
irregular NGDB sampling.

A further independently seeded study used irregular sites and omitted 80 of
144 candidate sites. Confirmation coverage was 96% under independent uniform
selection and 87% with spatial correlation. Under a specified mechanism favoring
high-concentration sites, the mean estimate was about 75.7 for a target of 50,
and target coverage was 0% in 100 confirmation datasets. This is a synthetic
selection mechanism, not a quantified estimate of NGDB bias. It demonstrates
that declustering observed records cannot reconstruct omitted low-concentration
sites. Defining a target population and investigating the sampling design or an
explicit spatial/selection model are necessary before routine area inference.

Separate research scripts now evaluate score-based spatial calibration with
declustering weights, including surveys whose locations and weights change.
In synthetic tests the weighted linear-Gaussian calculation met the 95% coverage
target within its assumed range of spatial dependence. The censored transfer
failed under heavy censoring: confirmation
coverage was about 86–94%, with downward point-estimate bias. Moderate-censoring
coverage was about 96–99%. These are restricted simulation results, not validated
USGS confidence intervals. The calibration routines remain outside the public
API, and weighted routine confidence limits remain missing.

## Detection-limit (censoring) conventions

`parse_usgs_q()` interprets the USGS NGDB encoding:

| Input | Meaning |
|-------|---------|
| `qvalue` < 0 | left-censored (non-detect) at LOD = `abs(qvalue)` |
| `qvalue` ending `.1111` or `.11111` | right-censored at the upper limit of determination |
| `qualifier` `L`, `N`, or `<` | left-censored |
| `qualifier` `G` or `>` | right-censored |
| `qvalue` > 0 | detected value |

Some NGDB CSV exports round `qvalue` to two decimal places, hiding the
upper-censoring suffix. When the first best-ranked value in `summary` contains
an explicit suffix and rounds to the reported qvalue, the builders recover
that upper limit only when the first entry is well formed, a nonblank/non-UNKNOWN
method matches, and both encodings are usable and compatible. The original
`qvalue_raw` is preserved and `censoring_source` records the recovery.
`summary_recovery_status` in long builder output records recovery or its blocking
reason. Blocked recovery leaves the original qvalue/qualifier coding intact;
it does not borrow another method's limit or certify the original measurement.
This guard is new in1.2.0.9003 and can change previously unguarded recovery.
Ordinary values ending in `.11` are not inferred to be censored from qvalue alone.

The MLE then represents left-censored observations as the interval `(0, LOD]`,
right-censored as `(ULOD, Inf)`, and detections as exact, fitting a lognormal by
`survival::survreg(..., dist = "lognormal")`; the geometric mean is
`exp(meanlog)`.

## Function reference

**Parsing & building**
- `parse_usgs_q()` — USGS qvalue/qualifier → censoring fields
- `build_ag_from_chem()` / `build_ag_from_csv()` / `build_element()` — censoring-aware element tables
- `build_metals_wide()` / `build_ag_metals_state()` — wide metals tables (substitution)
- `read_clean_main()` — read & clean NGDB sample metadata
- `join_soil_main()` / `join_sed_main()` — join chemistry to metadata by `lab_id`
- `norm_state()` — normalise state names to USPS codes

**Estimation**
- `fit_cen_lognorm()` — left+right-censored lognormal MLE (geomean + CI)
- `mle_by_state()` — censoring-aware lognormal geomeans by state
- `compare_censored_distributions()` — AIC comparison of censored parametric distributions
- `spatial_coverage_by_state()` — coordinate counts and occupied grid cells by state
- `cell_declustering_weights()` — weights sharing influence within occupied projected cells
- `declustering_sensitivity()` — weighted censored estimates across cell sizes
- `add_substitution()` — LOD/2 and upper-limit substitution values
- `fit_lognorm_with_quants()` — censoring-aware lognormal MLE with fitted quantiles
- `robust_cor_table()` — bootstrap Spearman correlations of metals vs a target
- `robust_corr_vs_Ag()` — bootstrap Spearman correlations of metals vs silver (`Ag_mgkg`); used for the element–Ag correlation tables
- `run_pca_robust()` — robust log10 PCA on metals

**Visualisation**
- `plot_lognormal_fit()` — fitted censored lognormal density with geomean and CI
- `plot_ranked_bar()` — ranked barplot of state geomeans (log scale)
- `plot_us_map()` / `state_hotspot_map()` — choropleth / hotspot maps

## Data

The geoCensor package **does not include geochemical data**. Instead, it is designed to work with user-downloaded USGS National Geochemical Database exports, including the soil and stream-sediment “best value” and “chemistry” tables, together with the “main” sample-metadata table. These files are available through the [USGS Mineral Resources Online Spatial Data portal](https://mrdata.usgs.gov/ngdb/).


## Citation

If you use geoCensor, please cite it. The repository includes citation metadata in `CITATION.cff`:

> Edwards, R. (2026). *geoCensor: Censoring-Aware Analysis of USGS National
> Geochemical Database Data* (version 1.2.0). https://github.com/diskoisland/geoCensor

## License

MIT © 2026 Ross Edwards
