# geoCensor 1.1.0 release candidate

Prepared November 1, 2026. Version metadata is 1.1.0; no release tag is published.

- Scope: USGS censoring/data-integrity improvements, multi-element workflows,
  descriptive declustering/support diagnostics and empirical censoring bounds.
- Ordinary model intervals assume independent observations and an appropriate
  lognormal model. Weighted routine limits remain NA. Spatial inference is
  experimental and research routines remain outside the installed API.
- Compatibility changes: duplicate sample/element pairs must be selected before
  widening; quantile summaries now distinguish insufficient data and failed fits.
- The clean 1.1.0 source package passed local `R CMD build --no-manual` and
  `R CMD check --no-manual` (Status: OK), including vignette creation/rebuilding.
  The full suite passed 77 tests / 410 expectations, with no failures, errors,
  warnings or skips. Installed single/all-element command-line examples passed
  on synthetic data, retaining leading-zero IDs, unsupported-species statuses
  and missing weighted confidence limits.
- Cross-platform GitHub Actions checks must pass on the exact review-branch
  commit before tagging. No merge, release or tag is performed by preparation.
- Local output artifacts and source data are not included in the publication
  allowlist. Build-log output links refer to unbundled local research artifacts.

The workflow follows the [r-lib standard check template](https://github.com/r-lib/actions/blob/v2/examples/check-standard.yaml).

Checked source archive SHA-256:

`9c696f14e9951160e9e1921afdb6897e22ce5531dcd809f7c5b73da4934b2767`

The checked archive has 101 files, generated from 97 explicitly selected source
files. The repository publication allowlist also includes CI, citation metadata
and development research sources. No tag or release is created by this step.
