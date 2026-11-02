# geoCensor 1.1.0 release validation

Prepared November 1, 2026. Release record and downloadable assets:
[v1.1.0](https://github.com/diskoisland/geoCensor/releases/tag/v1.1.0).

## Scope and compatibility

USGS censoring/data-integrity improvements, multi-element workflows,
descriptive declustering/support diagnostics and empirical censoring bounds.
Duplicate sample/element pairs must be selected before widening; quantile
summaries now distinguish insufficient data and failed fits. Default successful
quantile output names/order and concentration units are preserved.

Ordinary model intervals assume independent observations and an appropriate
lognormal model. Weighted routine limits remain NA. Spatial inference is
experimental, and research routines remain outside the installed API.

## Validation record

The candidate passed 77 tests / 410 expectations with no failures, errors,
warnings or skips. Installed single/all-element examples passed on synthetic
data, retaining leading-zero IDs, unsupported/all-left statuses and missing
weighted confidence limits. The final release archive passed a fresh
`R CMD build --no-manual` and `R CMD check --no-manual` (Status: OK), including
vignette creation/rebuilding and the full package test suite.

The final package's R, help, tests, installed examples and vignette sources are
identical to the checked candidate. Final installed-package changes are only
the README/NEWS release wording; citation metadata adds the release date.
The final source stage contains 97 files and its archive 101 files, including
generated vignette material. The archive excludes PFAS, source data, papers,
research scripts and generated research outputs.

Archive SHA-256:

`bfd47d1d89828b8cd87efe04a835cdbddb241cefaa12112d9e919fb1cb1232e2`

## Review, checks and publication gates

[PR #1](https://github.com/diskoisland/geoCensor/pull/1) records the review.
Both [candidate](https://github.com/diskoisland/geoCensor/actions/runs/41385826608)
and [documentation-head](https://github.com/diskoisland/geoCensor/actions/runs/41388359960)
matrices passed all five jobs: macOS/current R, Windows/current R,
Ubuntu/current R, Ubuntu/oldrel-1 and Ubuntu/devel.

Finalization requires the final metadata head's five checks to pass, then a
normal merge matching that exact PR head. The merged tree must match the checked
release tree and all five main-branch checks must pass before tagging.
Publication uses an annotated `v1.1.0` tag pinned to that merge commit, with no
force update or deletion. The release includes the checked archive, checksums
and an exact-commit validation receipt. Downloaded asset bytes are verified.

The release's `VALIDATION.md` asset records the final head, merge commit,
workflow links and archive checksum. Local research artifacts remain unbundled;
links into `output/` in the development log refer to those local records.


## Workflow input correction

Before tagging, the initial main run exposed an unsupported `check_args`
workflow input. The official action takes `args`; its default already matched
our intended `c("--no-manual", "--as-cran")`, so earlier checks ran as intended.
The input name was corrected through a separate review change. Final publication
is gated on green checks of the corrected head and its merged main commit.
The checked package files/archive are unchanged by this configuration-only fix.
