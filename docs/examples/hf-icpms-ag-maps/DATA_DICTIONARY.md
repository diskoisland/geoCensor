# Cohort structure and provenance

`data.RDS` is a base-R list. It contains data and metadata, with no executable
functions, spatial boundary objects, local filesystem paths or private campaign
records. Missing source values remain `NA`.

| Component | Content |
| --- | --- |
| `chemistry` | All 119,490 selected original Ag chemistry reports and their separately recorded interpretations |
| `metadata` | Unique referenced geology/sample rows, selected fields only |
| `methods` | Source method descriptions, digestion descriptions and method-reference keys |
| `snapshot` | Source label, table basenames, row counts, byte sizes and SHA256 hashes |
| `selection` | Method/media eligibility and explicit absence of nominal-limit or preparation filters |
| `interpretation` | Qualifier, unit, role, date and attribution qualifications |
| `citation`, `rights` | Source citation and data licensing |

## Original chemistry and interpreted observations

| Field | Meaning |
| --- | --- |
| `national_chem_row_ordinal_1based` | One-based row in the original national chemistry table; the report identity |
| `metadata_geol_row_ordinal_1based` | Exact linked row in the original national geology table |
| `lab_id`, `group_id` | Source identifiers; not guaranteed independent sampling events |
| `parameter`, `species`, `units` | Original source labels; species Ag, units ppm |
| `analytic_method`, `lab_name` | Original method and laboratory codes |
| `qualified_value`, `qualifier` | Unchanged original numeric encoding and qualifier |
| `parse_status`, `parse_rule` | Reviewed interpretation recorded separately from source fields |
| `result_value` | Positive numerical-coded result, present only for numerical observations |
| `reporting_lower_bound` | Actual lower bound for a right-censored observation |
| `reporting_upper_bound` | Actual upper bound for a left-censored observation |
| `medium_class`, `state_metadata` | Recorded soil/sediment class and linked jurisdiction |
| `metadata_join_status`, `jurisdiction_class`, `sediment_subtype` | Explicit linkage and source-derived classifications |

Numerical results are not automatically certified as accurate. `N` means not
detected; `L` means detected below the lower determination limit; `<` indicates
a below-limit result without specifying detection state. They are retained
separately although all enter this model as left-censored observations. `>`
and `G` identify above-limit results. Actual bounds are preserved from the
reviewed source encoding, including its special upper-limit suffix forms.
Nominal method limits are not used to assign a new bound or recensor a result.

## Sample metadata

The link uses original geology-row ordinal, `lab_id` and `group_id`, with exact
agreement required. The retained fields include source `OBJECTID`, NURE and
field IDs, job and project labels, collection/submission dates, country/state,
quadrangle, district/deposit/mine names, latitude/longitude, medium/subtype,
collection method, depth/horizon/zone, organics, preparation and mesh-size text.
No preparation class is inferred from an empty field. `date_submitted` is not
an analysis date. Unlabelled `reanalysis` values are not renamed “original”.

Some chemistry reports refer to the same metadata row. The metadata table
stores each referenced source row once; chemistry reports are not collapsed,
averaged or silently removed. A metadata row, group ID or repeated coordinate
does not certify an independent site or duplicate field sample.

## Map and reference results

`reference-cells.RDS` stores separately verified `pooled`, `soil` and `sediment`
cell tables. Grid indices and projected cell midpoints identify cells. Counts
include reports by censoring direction, numerical values and coordinate pairs.
`GM_mgkg` is missing for a cell without a displayed estimate. The remaining
fields retain the conditional lognormal parameters and numerical diagnostics.
These tables contain derived summaries, not replacements for source observations.

`hf_observations()` joins the two cohort tables and selects numerical values or
actual bounds. It adds a CONUS map flag from recorded state and valid WGS84
coordinates. Alaska, Hawaii and territories remain in the full RDS but are not
included in these contiguous-U.S. figures. Invalid/missing coordinates do not
become (0,0). `hf_project()` creates new projected coordinates and grid indices;
original coordinates are unchanged.

The cohort can be traced back to the source snapshot using its table hashes
and row identifiers. This is reversibility of the extraction and interpretation,
not a claim that the original archive can be reconstructed from this subset.
