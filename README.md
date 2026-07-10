# geoCensor

**Censoring-aware analysis of USGS National Geochemical Database (NGDB) data.**

I wrote geoCensor for working with USGS National Geochemical Database tables where a lot of the useful information is tucked inside detection limits, qualifiers, and non-detects. Those details matter: in the NGDB, detection limits vary by method and over time, so dropping non-detects or replacing everything with LOD/2 can move the center of the data in ways that are hard to see later.

The package turns USGS `qvalue` and `qualifier` fields into concentration intervals, including open-ended intervals for non-detects and greater-than values. It then estimates geometric means and confidence intervals with maximum-likelihood lognormal models using interval-censored survival regression. It also includes LOD/2-style substitution summaries for comparison, distribution diagnostics, spatial joins, robust correlation analysis, robust PCA, and mapping helpers.

To generate the MLE:
1. non-detects were coded as `(0, LOD]`;
2. greater-than values were coded as `(ULOD, Inf)`; and
3. measured detections were treated as exact concentrations.

The lognormal model is fit with `survival::survreg(..., dist = "lognormal")`, and `exp(meanlog)` is used as the geometric mean.


## Installation

```r
# install.packages("remotes")
remotes::install_github("diskoisland/geoCensor")
```

## Quick start

```r
library(geoCensor)

# 1. Parse a chemistry table with lab_id, species, unit, qvalue,
#    and qualifier columns into censoring fields for one element,
#    converting supported units to ug/kg
ag <- build_ag_from_csv("ngdbsoil/bestvalue.csv", element = "Ag", to_ugkg = TRUE)

# 2. Attach sample metadata (state, coordinates) by lab_id
main <- read_clean_main("ngdbsoil/main.csv")
ag_main <- join_soil_main(ag, main)

# 3. Censoring-aware lognormal MLE geometric mean, by state
gm_by_state <- mle_by_state(ag_main)
head(gm_by_state)        # state, n_total, n_left, n_right, geomean, ci_low, ci_high

# 4. Spatial coverage diagnostics for the same state metrics
coverage_by_state <- spatial_coverage_by_state(ag_main)
gm_by_state <- dplyr::left_join(gm_by_state, coverage_by_state, by = "state")

# 5. National pooled estimate
fit_cen_lognorm(ag_main)

# 6. Inspect the fitted lognormal model
plot_lognormal_fit(ag_main)

# 7. Map / rank
plot_ranked_bar(gm_by_state)
plot_us_map(gm_by_state)
```

## Example scripts

Installed example scripts show how to build Ag MLE state tables and hotspot
maps from NGDB `bestvalue.csv` files:

```r
system.file("examples", "build_ag_mle_state_df.R", package = "geoCensor")
system.file("examples", "ag_seed_state_mle.R", package = "geoCensor")
system.file("examples", "test_ag_spatial_coverage_step_by_step.R", package = "geoCensor")
system.file("examples", "ag_state_mle_spatial_coverage.R", package = "geoCensor")
system.file("examples", "test_ag_lognormal_fit_plot_step_by_step.R", package = "geoCensor")
system.file("examples", "ag_lognormal_fit_plot.R", package = "geoCensor")
system.file("examples", "ag_hotspot_map.R", package = "geoCensor")
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
- `compare_censored_distributions()` — AIC comparison of censored parametric distributions
- `spatial_coverage_by_state()` — coordinate counts and occupied grid cells by state
- `add_substitution()` — LOD/2 and upper-limit substitution values
- `fit_lognorm_with_quants()` — censoring-aware lognormal MLE with fitted quantiles
- `robust_cor_table()` — bootstrap Spearman correlations of metals vs a target
- `robust_corr_vs_Ag()` — bootstrap Spearman correlations of metals vs silver (`Ag_mgkg`); used for the element–Ag correlation tables
- `run_pca_robust()` — robust log10 PCA on metals

**Visualisation**
- `plot_lognormal_fit()` — fitted censored lognormal density with geomean and CI
- `plot_ranked_bar()` — ranked barplot of state geomeans (log scale)
- `plot_us_map()` / `state_hotspot_map()` — choropleth / hotspot maps

## Data

The geoCensor package **does not include geochemical data**. Instead, it is designed to work with user-downloaded USGS National Geochemical Database exports, including the soil and stream-sediment “best value” and “chemistry” tables, together with the “main” sample-metadata table. These files are available through the [USGS Mineral Resources Online Spatial Data portal](https://mrdata.usgs.gov/ngdb/).


## Citation

If you use geoCensor, please cite it. The repository includes citation metadata in `CITATION.cff`:

> Edwards, R. (2026). *geoCensor: Censoring-Aware Analysis of USGS National
> Geochemical Database Data* (version 0.1.0). https://github.com/diskoisland/geoCensor

## License

MIT © 2026 Ross Edwards
