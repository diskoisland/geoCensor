# geoCensor 1.2.0 release checks

Prepared November 4, 2026. Source version: **1.2.0**. This file describes the
release procedure; it is not itself evidence that a tag or release exists.
Published artifacts are listed on the [releases page](https://github.com/diskoisland/geoCensor/releases).
An earlier installed local1.2.0 candidate is distinct from the final release.

## Reviewed features

The candidate includes method/source-summary review, guarded upper-limit recovery,
sampling-support audits and position memberships, lossless chemistry import,
species-selective batch reading, disk export with provenance and output reservations,
and an Upper Illinois River Basin worked example. NEWS.md records the compatibility
changes. Original data, dates, methods, QA and dilution are not inferred.

Installed commands check final1.2.0 or newer and required exported capabilities.
They must accept an actual final-version installation and reject incomplete older
local candidates even if their version label is1.2.0. Restart R after upgrading.
Use the source.tar.gz installer; a review-bundle ZIP is not an R package installer.

## Illinois example

The vignette and installed illinois_case_study.R select the named project from
user-downloaded sediment main/bestvalue exports without chemistry filtering.
They preserve source text and separate record ordinals, append explicitly
projected coordinates, retain all reported species/unavailable fits, and save
support/results, figures and provenance. The optional sf package is needed only
for coordinate preparation; users may supply reviewed EPSG:5070 coordinates.
Source data and the private GIS/research collection are not bundled.

Before publication, execute both projection paths on the frozen real cohort in
an isolated final-version library. Verify439 metadata records,20,422 selected
chemistry records,48 species and332 recorded geographic positions for that exact
snapshot. Compare direct API and installed-command results, preserve all-left
species and missing weighted limits, verify source fields/indices and inspect
figures. Counts do not certify physical sites, specimens or independence.

## Exact-artifact validation

Build from an explicit source allowlist, check/rebuild both vignettes, run package
tests and execute installed commands using the resulting archive. Verify source
and archive fingerprints and obtain all five final-head CI jobs: Windows/macOS/
Linux release, Linux oldrel-1 and devel. Earlier checked development heads are
baseline evidence, not a receipt for release edits. Retain failures and amendments.

Audit public Git and package archives to exclude private development/literature
logs, fetch lists, Papers, raw data, generated research bundles, unrelated files
and local root Git history. The underlying review stack is PR4 -> PR6 -> PR7 ->
PR8 -> PR9 -> PR10; PR5's guard is already included in PR4. Verify live base/head
ancestry and reviews before integration. A release-integration branch must contain
the reviewed final tree without unreviewed scope or data changes.

## Scientific limits

These software/numerical checks do not validate geological background, anomaly
classes, ore probability, human/geological attribution or spatial confidence
coverage. Ordinary intervals assume independent observations and an appropriate
working model. Routine declustering-weighted limits remain unavailable and spatial
procedures remain experimental. Uncensored known/fitted covariance controls are
private research, not production confidence routines. Further covariance-refitting,
censoring and model-misspecification validation precedes any demonstration paper
promoting spatial inference. Current hydrology/context does not establish historical
sample provenance; missing information remains missing.

## Publication

Align version, citation, installation instructions and NEWS. Integrate the reviewed
stack and final release changes, verify that main has the reviewed tree and all
five checks pass on its exact commit, then pin an annotatedv1.2.0 tag to that
commit. Publish the checked source archive, concise release notes, checksums and
reproduction instructions. Download the release asset and verify its bytes.
Never move an existing tag or equate a local candidate with a published release.
