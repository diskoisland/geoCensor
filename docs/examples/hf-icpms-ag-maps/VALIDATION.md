# Reproduction checks

The fresh R replay on 12 November 2026 reproduced all **5,495 displayed cell
means**, with the same occupied and unavailable cells as the earlier,
independently checked calculations.

| View | Displayed means checked | Maximum relative GM difference |
| --- | ---: | ---: |
| Pooled soil and sediment | 2,616 | 6.20 × 10⁻⁸ |
| Soil | 974 | 3.98 × 10⁻⁸ |
| Sediment | 1,905 | 3.24 × 10⁻⁸ |

The prescribed tolerance is 1 × 10⁻⁴. The
[machine-readable receipt](validation.json) records the actual input hashes,
software versions and comparison results. All six final map images were
visually inspected for labels, shared scales, hollow cells and clipping.

The RDS preserves the selected original chemistry rows exactly. Checks cover
raw qualifiers and actual bounds, exact source-row/lab/group joins, method and
medium counts, all cell counts, projected geometry, likelihood derivatives,
uncensored and censored synthetic fits, and all-censored, constant-value and
insufficient-support cases. Values outside the display range retain endpoint
colours. No spatial confidence intervals are calculated.

Three pooled cells and one soil cell had no result from the released
geoCensor API's default numerical start. The separately recorded multistart
refinement recovered their conditional maximum-likelihood fits and matched
the archived results. No sediment cell needed that recovery. This numerical
recovery does not establish scientific reliability or resolve weak information
below reporting limits.

These checks establish software and calculation reproducibility for this
specific cohort. They do not validate regional representativeness, harmonize
sampling/preparation methods, demonstrate a lognormal model everywhere, or
establish sources of observed silver. The [validation guide](../../validation.md)
describes the wider scientific limits.
