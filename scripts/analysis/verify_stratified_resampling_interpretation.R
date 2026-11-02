#!/usr/bin/env Rscript
# Original deterministic interpretation checks, not author-code replication.
out <- "output/geochemical_database_papers_2026-11-01"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Uniformly selecting one observation in each occupied cell gives the same
# expected mean and empirical CDF as normalized occupied-cell weights.
a <- c(1, 5, 10, 100); b <- c(0, 4, 7)
draws <- as.matrix(expand.grid(a = a, b = b))
z <- c(a, b); w <- c(rep(1 / (2 * length(a)), length(a)),
                     rep(1 / (2 * length(b)), length(b)))
threshold <- sort(unique(z))
expected_cdf <- vapply(threshold, function(t) mean(rowMeans(draws <= t)), 0.0)
weighted_cdf <- vapply(threshold, function(t) sum(w[z <= t]), 0.0)
stopifnot(max(abs(expected_cdf - weighted_cdf)) < 1e-14,
          abs(mean(rowMeans(draws)) - sum(w * z)) < 1e-14)
write.csv(data.frame(threshold = threshold, expected_resample_cdf = expected_cdf,
  weighted_cdf = weighted_cdf), file.path(out, "resampling_weight_equivalence.csv"), row.names = FALSE)

# Each column is a valid one-from-A, one-from-B draw. Rankwise medians need
# not themselves be a realization containing one observation from each cell.
resamples <- rbind(A = c(1, 10, 5), B = c(4, 7, 0))
ranked <- apply(resamples, 2L, sort)
medians <- apply(ranked, 1L, median)
source_values <- c(1, 10, 5, 4, 7, 0)
source_cells <- c("A", "A", "A", "B", "B", "B")
median_cells <- source_cells[match(medians, source_values)]
stopifnot(identical(unname(medians), c(1, 5)), all(median_cells == "A"))
with_extremes <- unique(c(medians, min(source_values), max(source_values)))
extreme_cells <- source_cells[match(with_extremes, source_values)]
stopifnot(sum(extreme_cells == "A") == 3L, sum(extreme_cells == "B") == 1L)
write.csv(data.frame(value = with_extremes, cell = extreme_cells,
  in_rank_medians = with_extremes %in% medians),
  file.path(out, "rank_median_geometry_example.csv"), row.names = FALSE)
# The effect is not restricted to a particular three-draw random sequence:
# enumerating every equally probable pair also permits both medians from A.
pairs <- as.matrix(expand.grid(A = c(1, 9, 10), B = c(0, 4, 7)))
ranked_pairs <- t(apply(pairs, 1L, sort))
pair_medians <- apply(ranked_pairs, 2L, median)
stopifnot(identical(unname(pair_medians), c(1, 9)))
pair_values <- c(1, 9, 10, 0, 4, 7)
pair_cells <- c("A", "A", "A", "B", "B", "B")
appended <- unique(c(pair_medians, min(pair_values), max(pair_values)))
cell <- pair_cells[match(appended, pair_values)]
stopifnot(sum(cell == "A") == 3L, sum(cell == "B") == 1L)
write.csv(data.frame(value = appended, cell = cell,
  in_rank_medians = appended %in% pair_medians),
  file.path(out, "complete_pair_rank_medians.csv"), row.names = FALSE)
writeLines(c("Deterministic interpretation checks — November 1, 2026",
  "Exact enumeration verified the mean/CDF equivalence for uniform one-per-occupied-cell draws and normalized cell weights.",
  "Three valid two-cell draws yield rank medians 1 and 5, both from cell A.",
  "Appending missing extrema gives three observations from A and one from B; equal cell representation is not preserved.",
  "Complete enumeration of nine equally likely A={1,9,10}/B={0,4,7} pairs also yields two A rank medians and A3/B1 after extrema.",
  "This checks bookkeeping/geometry properties only, not the author's software or the empirical effectiveness of the method.",
  "Neither occupied-cell scheme supplies values for empty cells, independent new observations or confidence intervals."),
  file.path(out, "interpretation_checks.txt"))
invisible(file.copy("scripts/analysis/verify_stratified_resampling_interpretation.R", out, overwrite = TRUE))
cat("Exact resampling/weight equivalence and rank-median geometry checks passed.\n")
