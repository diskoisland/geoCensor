# Development analyses

These scripts retain geoCensor's USGS audits and statistical-method research.
They are excluded from the installed R package and are not production APIs.

Use the installed examples in `inst/examples` and the USGS workflow vignette
for ordinary package use. Research scripts can require additional R packages,
local NGDB files and saved simulation inputs that are not bundled here. Their
fixed input/output paths and study-specific fingerprints are intentional;
consult each study's protocol and source snapshots before attempting a rerun.
A completed study may refuse reruns if current sources differ from its frozen
signature. Preserve its recorded sources rather than replacing its results.

Detailed study records and `output/` artifacts are maintained locally and are
not bundled with this repository. Numerical convergence, likelihood preference and
sample-count guards do not establish valid spatial inference. Routine weighted
confidence limits remain missing. Confidence limits for the spatial research
remain experimental, and no real-NGDB site-block test is automatically run.

The release package's tests use synthetic data; running these private analyses
is not required to install geoCensor or to run its package test suite.

## Method-preserving South Carolina validation

The following observational study preserves all nine NGDB best-value fields,
including `amethod`, and the source data-record ordinal. Soil and sediment are
processed separately, with recorded collection types retained as groups. The
pilot uses soil Au/Al labelled single/grab; expansion processes every reported
species/property in explicitly U.S.-identified South Carolina records. This is
a pipeline reconciliation, not validation of spatial confidence coverage.

Run from the repository root, using a **new** output directory:

```sh
python3 scripts/analysis/prepare_usgs_method_validation.py output/usgs_method_validation_2026-11-02
Rscript scripts/analysis/run_usgs_method_validation.R output/usgs_method_validation_2026-11-02 pilot soil
Rscript scripts/analysis/run_usgs_method_validation.R output/usgs_method_validation_2026-11-02 expansion soil
Rscript scripts/analysis/run_usgs_method_validation.R output/usgs_method_validation_2026-11-02 expansion sed
python3 scripts/analysis/verify_usgs_method_validation.py output/usgs_method_validation_2026-11-02
Rscript scripts/analysis/verify_usgs_method_validation.R output/usgs_method_validation_2026-11-02
```

Prerequisites are Python 3 (standard library), R with `data.table`, `sf` and the
geoCensor dependencies, the checked installed geoCensor 1.1.0.9003 library at
`output/batch_method_overlap_2026-11-02/geoCensor.Rcheck`, and its canonical
`package_source_manifest.csv`. Local source inputs are `../data/ngdbsoil/` and
`../data/ngdbsed/` (`main.csv`, `bestvalue.csv`), plus existing projected metadata
and CRS definitions in `output/usgs_all_elements_2026-10-31/{soil,sed}/`. These
prerequisites and study artifacts are not distributed. To conduct a new study
with different prerequisites, first copy and adapt the preparation/driver
sources and declare the new protocol; do not alter a completed study's sources
or replace its results.

Preparation freezes the protocol, inputs and source fingerprints before reading
chemistry values. The driver independently reconciles the saved projection,
preserves per-species checkpoints, and resumes only with matching signatures.
The Python verifier re-streams the original records and independently checks
parsing, counts and every method-pair spatial overlap. The R verifier checks a
prespecified six-element panel against direct survival fits, independent grid
weights and empirical CDF bounds, plus analytic uncensored means. Missing
results and insufficient-data groups remain explicit. Ordinary model intervals
are conditional on unvalidated reporting/model/independence assumptions;
weighted routine limits remain missing. Geographic overlap does not establish
paired samples, method agreement, common sampling support or physical sites.
