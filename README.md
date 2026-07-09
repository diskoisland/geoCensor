# geoCensor

**Censoring-aware analysis of USGS National Geochemical Database (NGDB) data.**

Geochemical concentration data are heavily **left-censored**: a large fraction of
samples fall below the analytical detection limit, and the limit itself varies by
method and era. Naïve approaches (dropping non-detects, or substituting LOD/2)
bias estimates of central tendency. **geoCensor** treats detection limits
explicitly: it parses the USGS NGDB `qvalue`/`qualifier` encodings into
interval-censored fields and fits **maximum-likelihood lognormal models** (left-
and right-censored, via interval-censored survival regression) to estimate
**geometric means with confidence intervals**, optionally by state. It also
provides substitution-based alternatives for comparison, spatial joins, robust
correlations, robust PCA, and mapping.

## Installation

```r
# install.packages("remotes")
remotes::install_github("diskoisland/geoCensor")
```

## Quick start

```r
library(geoCensor)

# 1. Parse a chemistry table (lab_id, species, unit, qvalue, qualifier)
#    into censoring fields for one element, converting ppm -> ug/kg
ag <- build_ag_from_csv("ngdbsoil/bestvalue.csv", element = "Ag", to_ugkg = TRUE)

# 2. Attach sample metadata (state, coordinates) by lab_id
main <- read_clean_main("ngdbsoil/main.csv")
ag_main <- join_soil_main(ag, main)

# 3. Censoring-aware lognormal MLE geometric mean, by state
gm_by_state <- mle_by_state(ag_main)
head(gm_by_state)        # state, n_total, n_left, n_right, geomean, ci_low, ci_high

# 4. National pooled estimate
fit_cen_lognorm(ag_main)

# 5. Map / rank
plot_ranked_bar(gm_by_state)
plot_us_map(gm_by_state)
```

## Detection-limit (censoring) conventions

`parse_usgs_q()` interprets the USGS NGDB encoding:

| Input | Meaning |
|-------|---------|
| `qvalue` < 0 | left-censored (non-detect) at LOD = `abs(qvalue)` |
| `qvalue` ending `.1111` | right-censored at the upper limit of determination |
| `qualifier` `L` or `<` | left-censored | 
| `qualifier` `G` or `>` | right-censored |
| `qvalue` > 0 | detected value |

The MLE then represents left-censored observations as the interval `(0, LOD]`,
right-censored as `(ULOD, Inf)`, and detections as exact, fitting a lognormal by
`survival::survreg(..., dist = "lognormal")`; the geometric mean is
`exp(meanlog)`.

## Function reference

**Parsing & building**
- `parse_usgs_q()` — USGS qvalue/qualifier → censoring fields
- `build_ag_from_chem()` / `build_ag_from_csv()` / `build_element()` — censoring-aware element tables
- `build_metals_wide()` / `build_ag_metals_state()` — wide metals tables (substitution)
- `read_clean_main()` — read & clean NGDB sample metadata
- `join_soil_main()` / `join_sed_main()` — join chemistry to metadata by `lab_id`
- `norm_state()` — normalise state names to USPS codes

**Estimation**
- `fit_cen_lognorm()` — left+right-censored lognormal MLE (geomean + CI)
- `mle_by_state()` — censoring-aware lognormal geomeans by state
- `add_substitution()` / `fit_lognorm_with_quants()` — substitution-based comparison
- `robust_cor_table()` — bootstrap Spearman correlations of metals vs a target
- `robust_corr_vs_Ag()` — bootstrap Spearman correlations of metals vs silver (`Ag_mgkg`); used for the element–Ag correlation tables
- `run_pca_robust()` — robust log10 PCA on metals

**Visualisation**
- `plot_ranked_bar()` — ranked barplot of state geomeans (log scale)
- `plot_us_map()` / `state_hotspot_map()` — choropleth / hotspot maps

## Data

geoCensor ships **no data**. It operates on USGS National Geochemical Database
exports (soil and stream-sediment "best value" and "chemistry" tables, plus the
"main" sample-metadata table), available from the USGS Mineral Resources Online
Spatial Data portal (<https://mrdata.usgs.gov/geochemistry/>).

## Citation

If you use geoCensor, please cite it (see `CITATION.cff`):

> Edwards, R. (2026). *geoCensor: Censoring-Aware Analysis of USGS National
> Geochemical Database Data* (version 0.1.0). https://github.com/diskoisland/geoCensor

## License

MIT © 2026 Ross Edwards
