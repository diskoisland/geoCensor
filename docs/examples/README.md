# Reproduce the USGS examples

| Example | What it demonstrates |
| --- | --- |
| [National and state silver](national-state-silver.md) | Soil and sediment analyzed separately; cell declustering, grid sensitivity and unavailable estimates. Includes public aggregate tables and a standalone reproduction script. |
| [Upper Illinois metals](illinois-metals.md) | Project-level sediment summaries with censoring, spatial weights and multiple elements. |
| [Wisconsin metals](wisconsin-metals.md) | Semiquantitative methods and limited measurement information in a sediment project. |

## Sediment project workflow

The Illinois and Wisconsin pages use the released geoCensor 1.2.0 analysis
functions. Their additional reproduction script lives in this repository's
`docs/examples/` directory; it is not an installed 1.2.0 command.

## Inputs and selection

Download native sediment `main.csv` and `bestvalue.csv` exports separately from
the [USGS NGDB portal](https://mrdata.usgs.gov/ngdb/), retaining their provenance.
No geochemical database or raw selected-record table is bundled here. The tested
snapshot is identified by:

| Input | SHA256 |
| --- | --- |
| `main.csv` | `e67d5bfbebf6ceb9502e9e421948d1711c47537b51414f98a362fa0ad7d8faa9` |
| `bestvalue.csv` | `291591c10a664b5be7fca534bdb746d01aa9306f3280f2e48ad36641cbeff60c` |

The script selects these two cohorts from metadata before examining chemistry:

| Output folder | Recorded state | Exact project label |
| --- | --- | --- |
| `illinois/` | `IL` | `NAWQA Upper Illinois R. Basin` |
| `wisconsin/` | `WI` | `Wisconsin Lead-Zinc` |

All chemistry for the selected IDs is retained. Selection does not depend on
concentration, detection, fitted estimates or appearance in a figure. Both
cohorts must be present, with unique nonblank metadata IDs and the recorded
sediment medium. Unexpected media or malformed CSV records stop the run.

Counts can differ with another source version. Hashes identify a file snapshot;
they do not certify data accuracy. Keep your source version and its identifiers.

## Run

Install [geoCensor 1.2.0](https://github.com/diskoisland/geoCensor/releases/tag/v1.2.0)
and restart R. This example additionally uses `sf` to transform the published
WGS84 coordinates to EPSG:5070 metres:

```r
install.packages("sf") # only if needed
```

From a downloaded or cloned repository, run:

```sh
mkdir -p output
Rscript docs/examples/sediment-metals.R /path/to/ngdbsed/main.csv /path/to/ngdbsed/bestvalue.csv output/sediment_examples
```

Use your own quoted paths where they contain spaces. The output directory must
not already exist; its parent must exist. Keep output under an ignored local
directory, because it contains selected original-source records. Another run
uses another directory name; partial outputs from a write failure are retained.
The large CSVs are scanned in batches, including unselected records. Retained
memory grows with the selected cohorts rather than the entire chemistry export.

## Outputs and checks

Each cohort gets original selected metadata and chemistry, a projected analysis
copy, source-record indices, all-species estimates, method and support audits,
CDF bounds, selected-element summaries and figures. A serialized result preserves
the analysis objects. The root output records input/output MD5 fingerprints,
selection/configuration, projection software and session information. Input
checksums are verified again before output creation. Source-record ordinals
exclude the header and skipped blank lines.

The displayed six elements are Pb, As, Cd, Cu, Hg and Zn. The raw analysis retains
all 48 Illinois/32 Wisconsin species, including unavailable fits. A species absent
from another snapshot has no fabricated estimate. Cell sizes are fixed at 5, 10
and 20 km with 100 deterministic grid placements. Weighted and unweighted fits
use the same usable chemistry/coordinate subset. Output estimates use ug/kg;
figures convert to mg/kg. No spatial confidence limits are computed.

`lead_declustering.png`, `reporting.png`, `metals_sensitivity.png` and
`reporting_bounds.png` are generated where the required observations exist.
The public Wisconsin page emphasizes reporting and method review. Its Pb model
figure remains in the user's output for inspection, with the semiquantitative
warning; it is not used for a state-to-state comparison.

These are traceable descriptive examples. The [validation guide](../validation.md)
explains the remaining modeling, measurement and spatial-inference limits.

[Illinois](illinois-metals.md) · [Wisconsin](wisconsin-metals.md) · [Project](../../README.md)
