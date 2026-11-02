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
