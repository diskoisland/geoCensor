# geoCensor 1.2.0 release candidate validation

Prepared November 2, 2026. Candidate review: [PR #4](https://github.com/diskoisland/geoCensor/pull/4).
The checked zero-scale guard from [PR #5](https://github.com/diskoisland/geoCensor/pull/5)
is combined into the candidate. This checklist describes preparation, not a
published release. Published assets are on the [releases page](https://github.com/diskoisland/geoCensor/releases).

## Scope and compatibility

Method-review propagation, optional within-method summaries, geographic method
coverage/overlap, safe output-directory handling, lossless large-grid keys and
explicit zero-scale fitting-boundary rejection. See NEWS.md for additive output
fields and interpretation. The installed commands require version 1.2.0.

Raw records and empirical bounds remain available when parametric fitting is
unavailable. No reporting-class precision or missing metadata is inferred.
Ordinary limits depend on unvalidated model/independence assumptions; weighted
routine limits remain NA and spatial limits remain experimental. Method-specific
declustering and other research routines remain outside the installed API.

## Evidence and final checks

The integrated guard passed 123 tests / 831 expectations, no failures/errors/
warnings/skips, plus a full local package check including vignettes and all five
platform/R jobs. Independent observational studies retained 95,323 SC results;
the guarded method-specific study reconciled 7,035 model cases and every weight.
These are numerical/source checks, not validation of statistical interval coverage.

The final release version, installed commands and polished vignette require a
fresh exact-source archive build/check, installed single/all-element smoke checks,
archive-content/source fingerprints and five final-head platform checks. Detailed
receipts and archive SHA256 are recorded with the prepared artifact outside the
package; no passing status is inferred from an older head. The archive excludes
private working notes, papers, raw data, generated outputs and research scripts.

## Publication gates

Before any tag, confirm the final reviewed head and exact archive match the
validation receipt. A normal merge to main must preserve the checked tree, and
all five checks on that main commit must pass. Only then create a new annotated
version tag pinned to that commit and publish the checked archive, checksums and
exact-commit validation receipt. Never force-update an existing version tag.
Verify downloaded release asset bytes. This preparation step makes no main merge,
tag or release. Private development/literature notes remain local and untracked.

## Previous release

The unchanged [v1.1.0 release](https://github.com/diskoisland/geoCensor/releases/tag/v1.1.0)
retains its original validation assets. Its source archive SHA256 is
bfd47d1d89828b8cd87efe04a835cdbddb241cefaa12112d9e919fb1cb1232e2.
