# geoCensor development candidate and publication checks

Reviewed November 3, 2026. Current source: **1.2.0.9008 development**. The latest
published release verified during this audit is
[v1.1.0](https://github.com/diskoisland/geoCensor/releases/tag/v1.1.0);
main is a separate ref and does not contain the stacked development candidate.
An installed local 1.2.0 candidate does not establish a public 1.2.0 release.

## Review chain

The existing open stack is [PR4](https://github.com/diskoisland/geoCensor/pull/4)
into main, [PR6](https://github.com/diskoisland/geoCensor/pull/6) into PR4,
[PR7](https://github.com/diskoisland/geoCensor/pull/7) into PR6,
[PR8](https://github.com/diskoisland/geoCensor/pull/8) into PR7, and
[PR9](https://github.com/diskoisland/geoCensor/pull/9) into PR8.
The zero-scale guard from PR5 is already included in PR4. PR9's feature head is
`df40fcc283ef2554d2609f28c9c6afdff5087a2c`; it contains the earlier checked changes.
Resolve reviews in dependency order and verify base/head ancestry before any
merge. A newer documentation revision needs its own exact-source check receipt.
This checklist is preparation, not a merge or release record.

## Scope and installed workflows

The candidate includes method/source-summary review, guarded upper-limit recovery,
sampling-support audits and position memberships, lossless chemistry import,
optional species-selective batch reading and optional disk export with provenance
and owned output reservations. Existing analysis defaults are preserved. NEWS.md
records compatibility changes, including recovery guards that can change earlier
unguarded interpretations.

Minimum installed versions: sampling audit 1.2.0.9005; all-element command
1.2.0.9006; single-element command 1.2.0.9007; disk export command 1.2.0.9008.
Use the explicit development revision or checked source `.tar.gz` in README.md.
Restart R after upgrading a loaded package and compare installed/loaded versions.
The review ZIP is evidence, not an installer. Citation metadata identifies the
development source, while the published v1.1.0 citation remains on its tag.

## Validation evidence and limits

PR9 passed local source-package checking and all five exact-head CI jobs with
1,262 expectations and no failures, warnings, skips or notes. The full native
sediment export independently matches 12,098,896 records, 81 species labels,
108,890,064 original fields and every source ordinal. All 12 components of the
48-species pilot analysis are unchanged after export/reimport. Counts are chemistry
records and recorded positions, not certified physical sites or independent data.

The vignette contains a runnable synthetic export/read/analyze example, including
missing and unavailable results. For a final reviewed revision, build a fresh
allowlisted source archive, check/rebuild vignettes, execute installed workflows,
verify archive/source fingerprints and confirm all five final-head CI logs.
Keep receipts and checksums with the checked artifact; do not borrow passing
status from an earlier head. Exclude private notes, papers, raw data, generated
research outputs and root local history from public source and package archives.

These are software and numerical checks. Ordinary confidence limits still depend
on model/independence assumptions; weighted routine limits remain missing and
spatial procedures remain experimental. Unknown dates, sample support, methods,
QA, dilution and unsampled terrain remain unknown. Spatial-calibration research
stays outside the production API pending its own validation.

## Publication sequence

Agree on the final release version and reviewed tree, update version/citation/
installation metadata together, and validate its exact archive. A main merge must
preserve the reviewed tree, with all five checks passing on the resulting main
commit. Only then create a new annotated tag pinned to that commit and publish
the checked archive, checksums and exact-commit receipt. Never move an existing
version tag. Verify downloaded release-asset bytes. No main merge, tag or release
is performed by this readiness pass.
