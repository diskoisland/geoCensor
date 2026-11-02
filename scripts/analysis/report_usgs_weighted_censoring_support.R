#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
out <- "output/usgs_weighted_censoring_support_2026-11-01"
if (!file.exists(file.path(out, "COMPLETED.txt"))) stop("Support audit incomplete")
x <- rbindlist(lapply(c("soil", "sed"), function(medium)
  fread(file.path(out, medium, "weighted_censoring_support.csv"))))
stopifnot(nrow(x) == 22431L, all(x$n_exact + x$n_left + x$n_right == x$n_total))
valid <- x$n_total > 0
stopifnot(all(is.na(x$weighted_exact_fraction[!valid])),
  all(abs(x$weighted_exact_fraction[valid] + x$weighted_left_fraction[valid] +
    x$weighted_right_fraction[valid] - 1) < 1e-12))
summary <- x[, .(species = uniqueN(element), rows = .N, usable_rows = sum(n_total > 0),
  empty_rows = sum(n_total == 0), raw_below_5_percent = sum(exact_fraction < .05, na.rm = TRUE),
  weighted_below_5_percent = sum(weighted_exact_fraction < .05, na.rm = TRUE),
  raw_above_but_weighted_below = sum(n_exact >= 10 & exact_fraction >= .05 &
    weighted_exact_fraction < .05, na.rm = TRUE)), by = .(medium, cell_size)]
fwrite(summary, file.path(out, "support_summary.csv"))
review <- x[n_exact >= 10 & exact_fraction >= .05 & weighted_exact_fraction < .05]
setorder(review, weighted_exact_fraction, medium, element, state, cell_size)
fwrite(review, file.path(out, "weighted_support_review_cases.csv"))
unique_review <- unique(review[, .(medium, element, state)])
element_summary <- x[n_total > 0, .(state_cell_rows = .N,
  min_raw_exact_fraction = min(exact_fraction),
  min_weighted_exact_fraction = min(weighted_exact_fraction),
  median_weighted_exact_fraction = median(weighted_exact_fraction),
  below_weight_review_convention = sum(weighted_exact_fraction < .05)), by = .(medium, element)]
fwrite(element_summary, file.path(out, "element_support_summary.csv"))

png(file.path(out, "exact_weight_fractions.png"), width = 1850, height = 1100, res = 150)
layout(matrix(c(1, 2, 3, 3), 2L, 2L, byrow = TRUE), heights = c(1, .12))
par(mar = c(5, 5, 3, 1))
colors <- c("#2166ac", "#b35806", "#4d9221")
for (medium in c("soil", "sed")) {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Exact fraction of usable records",
    ylab = "Exact fraction of declustering weight", main = medium)
  abline(0, 1, col = "grey65", lty = 3)
  abline(h = .05, v = .05, col = "grey75", lty = 2)
  for (i in seq_along(c(5000, 10000, 20000))) {
    matching <- which(x$medium == medium & x$cell_size == c(5000, 10000, 20000)[i] & x$n_total > 0)
    z <- x[matching]
    points(z$exact_fraction, z$weighted_exact_fraction, pch = 16, cex = .35,
      col = adjustcolor(colors[i], alpha.f = .3))
  }
}
par(mar = rep(0, 4)); plot.new()
legend("center", legend = c("5 km cells", "10 km cells", "20 km cells"),
  col = colors, pch = 16, horiz = TRUE, bty = "n")
invisible(dev.off())
table <- vapply(seq_len(nrow(summary)), function(i) {
  a <- summary[i, ]
  sprintf("| %s | %d | %d | %d | %d | %d |", a$medium, a$cell_size / 1000,
    a$usable_rows, a$raw_below_5_percent, a$weighted_below_5_percent,
    a$raw_above_but_weighted_below)
}, "")
writeLines(c("# USGS censoring-weight support audit — November 1, 2026", "",
  "Recomputed descriptive geometry/censoring measures for all 80 soil and 81 sediment species in the completed U.S. analysis. Chemistry and cached projected metadata were read only. Existing concentrations, intervals and prior outputs were not refitted or rewritten.", "",
  "Every one of the 22,431 element/state/cell rows reconciled usable/excluded counts, left/right counts and Kish weight n with the original 5/10/20 km, 100-origin analysis. A direct public-API spot check was also saved for each medium. Source/input hashes and checkpoints are retained.", "",
  "## Measures", "",
  "exact_fraction counts exact observations within the same usable chemistry/coordinate subset. weighted_exact_fraction is their share of total declustering weight; weighted_left_fraction and weighted_right_fraction complete the shares. These are geometry-weight fractions, not independent information counts, known inclusion probabilities or tests of reliability. Empty subsets have missing fractions.", "",
  "## Descriptive review comparison", "",
  "The 5% share and ten-exact minimum below reuse the existing count-review conventions for comparison. They are not statistically validated thresholds, and no estimate/observation is excluded. State/cell rows for the same element share observations and are not independent replicates.", "",
  "| Medium | Cell km | Usable rows | Raw exact <5% | Weighted exact <5% | Raw clears review but weighted <5% |",
  "| --- | ---: | ---: | ---: | ---: | ---: |", table, "",
  sprintf("There are %d cell/element/state rows where raw counts have at least ten exact observations and at least 5%% exact values, while exact values carry under 5%% of the declustering weight. They represent %d distinct medium/element/state combinations across %d species labels. The complete retained review list is weighted_support_review_cases.csv.",
    nrow(review), nrow(unique_review), uniqueN(review$element)), "",
  "![Raw versus weighted exact fractions](exact_weight_fractions.png)", "",
  "These comparisons reveal how observation support changes after geometric weighting. They do not measure actual NGDB selection bias, establish a population mean, select a preferred cell size or certify a fitted distribution. Heavy censoring, sparse exact support, sample-support mixtures and preferential selection remain separate limitations.", "",
  "## Reproduction", "",
  "Run scripts/analysis/audit_usgs_weighted_censoring_support.R, then report_usgs_weighted_censoring_support.R. The audit recalculates chemistry validity and weights without concentration-model refits. Its manifest rejects changed inputs or sources. Public declustering_sensitivity() and analyze_elements() now expose the same descriptive fractions; package tests and source-package validation are separate from this all-data audit."),
  file.path(out, "README.md"))
cat("All-data support report reconciled:", nrow(review), "review rows across",
  nrow(unique_review), "medium/element/state combinations.\n")
