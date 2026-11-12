# Data licensing and attribution

Checked 12 November 2026.

## Selected geochemical data

The [USGS release landing page](https://www.usgs.gov/data/geochemical-data-rock-sediment-soil-and-mineral-samples-united-states-and-territories-1962)
for [DOI 10.5066/P1FXE9GK](https://doi.org/10.5066/P1FXE9GK) explicitly marks the
work **CC0 1.0 Universal**. The
[CC0 terms](https://creativecommons.org/publicdomain/zero/1.0/) permit copying,
modification and redistribution. The selected records in `data.RDS` retain that
status; the MIT software license does not replace the data rights statement.
Credit the USGS and cite the release as given in the example README.

The [USGS copyright policy](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits)
also describes USGS-produced data as U.S. public-domain information, while
distinguishing third-party copyrighted material. This example distributes
selected published chemistry and metadata values, rather than copying report
PDFs, photographs, vendor documents, agency logos or private laboratory data.
No third-party copyright marking was identified on the selected archive fields
or the release's data-rights statement.

The actual extract is from the supplied November 2025 Parquet snapshot.
Source-table hashes and identifiers are retained; the current 2026 catalogue
metadata is not a claim that a later download has identical bytes. Reformatting
as RDS and recording interpreted bounds do not replace original source fields
or imply USGS endorsement of the analysis.

## Code and state boundaries

Original example R code and documentation use the repository's
[MIT license](../../../LICENSE.md). Numerical reference summaries are supplied
for reuse with the CC0 source data.

State boundaries are drawn at runtime using the
[`maps` package](https://cran.r-project.org/package=maps), a separately licensed
**GPL-2** dependency, with the version recorded in each run. Boundary data are
not in `data.RDS` or `reference-cells.RDS`; this example does not relicense them
as CC0. Preserve boundary attribution and applicable GPL terms when reusing the
map layer. The package source and license remain available through
[CRAN](https://cran.r-project.org/src/contrib/Archive/maps/).
The package's `state` documentation identifies the underlying boundary source
as U.S. Department of the Census data; the software package's GPL-2 status is
recorded separately from the USGS data rights.

The plotted Ag values, state outlines, R code and related source publications
have distinct provenance. Availability of one does not grant blanket rights
to illustrations or other content from the cited publications.
