#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); pkgload::load_all(".", quiet = TRUE)})
setDTthreads(2L)
out <- "output/censored_cdf_diagnostics_2026-11-01"
if (dir.exists(out) && length(list.files(out))) stop("Use a fresh demonstration directory")
dir.create(out, recursive = TRUE)
source <- file.path("..", "data", "ngdbsoil", "bestvalue.csv")
metadata_file <- "output/usgs_all_elements_2026-10-31/soil/metadata_projected.csv"
audit_file <- "output/usgs_weighted_censoring_support_2026-11-01/soil/weighted_censoring_support.csv"
sources <- c("scripts/analysis/demonstrate_usgs_cdf_bounds.R",
  sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)))
files <- c(sources, source, metadata_file, audit_file)
hash <- tools::md5sum(files)
write.csv(data.frame(file = files, md5 = unname(hash)), file.path(out, "source_input_manifest.csv"), row.names = FALSE)
dir.create(file.path(out, "source_snapshot"))
invisible(file.copy(sources, file.path(out, "source_snapshot")))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
metadata <- fread(metadata_file, colClasses = c(lab_id = "character"))
metadata <- metadata[state == "SC"]
chemistry <- fread(source, select = c("lab_id", "species", "qvalue", "unit", "summary"),
  colClasses = c(lab_id = "character", qvalue = "character"))
chemistry <- chemistry[species %in% c("Au", "Al") & lab_id %in% metadata$lab_id]
audit <- fread(audit_file)
plots <- list(); cdfs <- list(); supports <- list()
for (element in c("Au", "Al")) {
  z <- build_element(as.data.frame(chemistry), element)
  ids <- match(trimws(as.character(z$lab_id)), metadata$lab_id)
  stopifnot(!anyNA(ids), !anyDuplicated(z$lab_id))
  z$x_m <- metadata$x_m[ids]; z$y_m <- metadata$y_m[ids]
  intervals <- geoCensor:::.censored_intervals(z, "detected_val", "left_cen", "left_limit",
    "right_cen", "right_limit", "usgs_cdf_demonstration")
  usable <- intervals$ok & is.finite(z$x_m) & is.finite(z$y_m)
  d <- z[usable, , drop = FALSE]
  d$cell_weight <- cell_declustering_weights(d, "x_m", "y_m", 10000, 10L)
  values <- c(intervals$lower[usable], intervals$upper[usable])
  values <- values[is.finite(values) & values > 0]
  thresholds <- 10^seq(log10(min(values)) - .25, log10(max(values)) + .25, length.out = 180L)
  ref_rows <- which(audit$element == element & audit$state == "SC" & audit$cell_size == 10000)
  reference <- audit[ref_rows]
  stopifnot(nrow(reference) == 1L, nrow(d) == reference$n_total)
  for (weighted in c(FALSE, TRUE)) {
    key <- paste(element, if (weighted) "cell10km" else "equal", sep = "_")
    result <- censored_cdf_bounds(d, thresholds,
      weight_col = if (weighted) "cell_weight" else NULL)
    if (weighted) stopifnot(
      abs(result$support$weighted_exact_fraction - reference$weighted_exact_fraction) < 1e-12,
      abs(result$support$weighted_left_fraction - reference$weighted_left_fraction) < 1e-12,
      abs(result$support$weighted_right_fraction - reference$weighted_right_fraction) < 1e-12)
    cdfs[[key]] <- cbind(element = element, state = "SC", weighting = key, result$cdf)
    supports[[key]] <- cbind(element = element, state = "SC", weighting = key,
      source_n_input = nrow(z), source_n_excluded = sum(!usable), result$support)
    plots[[key]] <- plot_censored_cdf_bounds(result, title = paste0(element,
      ", South Carolina soil\n", if (weighted) "10 km cell weights" else "Equal record weights"))
    saveRDS(result, file.path(out, paste0(key, ".rds")))
  }
}
fwrite(rbindlist(cdfs), file.path(out, "empirical_cdf_bounds.csv"))
fwrite(rbindlist(supports), file.path(out, "observation_support.csv"))
png(file.path(out, "usgs_cdf_bounds.png"), width = 2400, height = 1800, res = 180)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2L, 2L)))
for (i in seq_along(plots)) print(plots[[i]],
  vp = grid::viewport(layout.pos.row = if (i <= 2L) 1L else 2L,
                      layout.pos.col = if (i %% 2L) 1L else 2L))
invisible(dev.off())
saveRDS(plots, file.path(out, "plots.rds"))
stopifnot(identical(tools::md5sum(files), hash))
writeLines(c("# Censoring-aware empirical CDF demonstration — November 1, 2026", "",
  "Calculated empirical CDF bounds for all usable South Carolina soil Au and Al records in the existing USGS data, separately with equal weights and 10 km cell weights (100 origins). The same usable chemistry/coordinate subset is used within each element. Projected metadata and chemistry are read only; prior results are unchanged.", "",
  "Au was selected because the previous support audit identified little exact-observation weight; Al provides a largely observed comparison. This is a deliberately selected diagnostic example, not an independent performance study. A state/medium subset can still mix geology, sample supports and analytical methods; no homogeneous geological population is asserted.", "",
  "Each vertical segment shows the minimum and maximum empirical cumulative fraction at one threshold. Every usable record contributes; censored concentrations are not replaced by limits or limit fractions. No interpolation between evaluation thresholds, distribution fit, p-value or confidence interval is calculated. Zero censoring width does not imply zero population or spatial uncertainty.", "",
  "The weighted censoring shares reconcile to the prior all-element audit. Source/input hashes and source snapshots are retained. observation_support.csv includes excluded counts and weighted shares; empirical_cdf_bounds.csv retains all evaluations.", "",
  "![Observed-data censoring bounds](usgs_cdf_bounds.png)", "",
  "Use the plot to distinguish unknown distribution shape caused by limits from observed shape. Model adequacy requires additional domain review and properly calibrated diagnostics; a fitted population CDF outside these sample bounds is not an automatic rejection. The recent scale/spatial validation remains a separate requirement."),
  file.path(out, "README.md"))
writeLines("Real-data diagnostic examples complete; no concentration estimates/intervals refitted.", file.path(out, "COMPLETED.txt"))
cat("South Carolina Au/Al empirical censoring-bound examples complete.\n")
