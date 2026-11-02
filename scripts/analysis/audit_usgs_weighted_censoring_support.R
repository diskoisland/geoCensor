#!/usr/bin/env Rscript
# Descriptive all-element audit. Recomputes geometry weights, not estimates.
suppressPackageStartupMessages({
  library(data.table)
  pkgload::load_all(".", quiet = TRUE)
})
setDTthreads(2L)
out <- "output/usgs_weighted_censoring_support_2026-11-01"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sources <- c("scripts/analysis/audit_usgs_weighted_censoring_support.R",
  sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)))
input_files <- unlist(lapply(c("soil", "sed"), function(medium) c(
  file.path("..", "data", paste0("ngdb", medium), "bestvalue.csv"),
  file.path("output/usgs_all_elements_2026-10-31", medium, "metadata_projected.csv"),
  file.path("output/usgs_all_elements_2026-10-31", medium, "state_declustering.csv"))))
signature <- list(sources = tools::md5sum(sources), inputs = tools::md5sum(input_files),
  cell_sizes = c(5000, 10000, 20000), n_offsets = 10L, R = R.version.string)
manifest <- file.path(out, "audit_manifest.rds")
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), signature)) stop("Audit inputs or sources changed")
} else {
  if (length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Unverified files in audit directory")
  dir.create(file.path(out, "source_snapshot"))
  invisible(file.copy(sources, file.path(out, "source_snapshot")))
  saveRDS(signature, manifest)
  write.csv(data.frame(file = c(sources, input_files),
    md5 = unname(c(signature$sources, signature$inputs))),
    file.path(out, "source_input_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(c("USGS weighted-censoring support audit — November 1, 2026",
    "All 80 soil and 81 sediment species; U.S. records identified by cached projected metadata.",
    "State/domain subsets, 5/10/20 km cells and 100 deterministic origins match the completed all-element analysis.",
    "Recompute chemistry validity and geometry weights; reconcile every group count and Kish weight n with the original table.",
    "No new concentration estimate, confidence interval, automatic exclusion or preferred cell size.",
    "Raw and weighted exact/left/right fractions are descriptive. A 5% comparison uses an existing review convention, not a validated threshold.",
    "Source data and all previous results are read only. Checkpoints resume only with matching source/input/configuration fingerprints."),
    file.path(out, "AUDIT_PROTOCOL.txt"))
}
if (file.exists(file.path(out, "COMPLETED.txt"))) {
  cat("Completed support audit matches provenance; no results changed.\n"); quit(status = 0L)
}
spot_checks <- list()
for (medium in c("soil", "sed")) {
  folder <- file.path(out, medium)
  dir.create(folder, showWarnings = FALSE)
  dir.create(file.path(folder, "checkpoints"), showWarnings = FALSE)
  base <- file.path("output/usgs_all_elements_2026-10-31", medium)
  metadata <- fread(file.path(base, "metadata_projected.csv"),
    colClasses = c(lab_id = "character"))
  if (anyDuplicated(metadata$lab_id)) stop("Repeated cached metadata IDs")
  original <- fread(file.path(base, "state_declustering.csv"))
  inventory <- fread(file.path(base, "input_inventory.csv"))
  chemistry <- fread(file.path("..", "data", paste0("ngdb", medium), "bestvalue.csv"),
    select = c("lab_id", "species", "unit", "qvalue", "summary"),
    colClasses = c(lab_id = "character", qvalue = "character"))
  stopifnot(nrow(chemistry) == inventory$input_chemistry_rows)
  chemistry <- chemistry[lab_id %in% metadata$lab_id]
  stopifnot(nrow(chemistry) == inventory$us_chemistry_rows)
  setkey(chemistry, species)
  elements <- sort(unique(chemistry$species))
  stopifnot(length(elements) == inventory$n_species)
  spot_done <- FALSE
  for (i in seq_along(elements)) {
    element <- elements[i]
    checkpoint <- file.path(folder, "checkpoints", sprintf("%03d.rds", i))
    if (file.exists(checkpoint)) {
      cached <- readRDS(checkpoint); stopifnot(identical(cached$element, element))
      if (!is.null(cached$spot)) {spot_checks[[medium]] <- cached$spot; spot_done <- TRUE}
      next
    }
    records <- as.data.frame(chemistry[.(element)])
    if (anyDuplicated(records$lab_id)) stop("Duplicate chemistry sample/element")
    df <- build_element(records, element)
    id <- match(trimws(as.character(df$lab_id)), metadata$lab_id)
    if (anyNA(id)) stop("Chemistry failed cached-metadata match")
    df$state <- metadata$state[id]; df$x_m <- metadata$x_m[id]; df$y_m <- metadata$y_m[id]
    groups <- split(seq_len(nrow(df)), df$state)
    rows <- list(); k <- 0L; spot <- NULL
    for (region in names(groups)) {
      d <- df[groups[[region]], , drop = FALSE]
      intervals <- geoCensor:::.censored_intervals(d, "detected_val", "left_cen", "left_limit",
        "right_cen", "right_limit", "weighted_censoring_support_audit")
      usable <- intervals$ok & is.finite(d$x_m) & is.finite(d$y_m)
      z <- d[usable, , drop = FALSE]
      exact <- intervals$lower[usable] == intervals$upper[usable]
      left <- intervals$lower[usable] == 0; right <- is.infinite(intervals$upper[usable])
      for (size in signature$cell_sizes) {
        w <- cell_declustering_weights(z, "x_m", "y_m", size, signature$n_offsets)
        n <- length(w)
        fraction <- if (n) c(sum(w[exact]), sum(w[left]), sum(w[right])) / sum(w) else rep(NA_real_, 3L)
        if (n) stopifnot(all(is.finite(fraction)), all(fraction >= 0),
          abs(sum(fraction) - 1) < 1e-12)
        # Calculate indices outside data.table's column-name evaluation scope.
        matching <- which(original$element == element & original$state == region & original$cell_size == size)
        old <- original[matching]
        stopifnot(nrow(old) == 1L, n == old$n_total,
          sum(left) == old$n_left, sum(right) == old$n_right,
          nrow(d) == old$n_input, sum(!usable) == old$n_excluded)
        kish <- if (n) 1 / sum(w^2) else 0
        stopifnot(abs(kish - old$weight_effective_n) <= 1e-8 * max(1, kish))
        k <- k + 1L
        rows[[k]] <- data.frame(medium = medium, element = element, state = region,
          cell_size = size, n_input = nrow(d), n_excluded = sum(!usable), n_total = n,
          n_exact = sum(exact), n_left = sum(left), n_right = sum(right),
          exact_fraction = if (n) mean(exact) else NA_real_,
          weighted_exact_fraction = fraction[1L], weighted_left_fraction = fraction[2L],
          weighted_right_fraction = fraction[3L], weight_effective_n = kish)
        if (!spot_done && n >= 6L && n <= 100L) {
          public <- declustering_sensitivity(d, size, "x_m", "y_m")
          stopifnot(abs(public$weighted_exact_fraction - fraction[1L]) < 1e-12,
            abs(public$weighted_left_fraction - fraction[2L]) < 1e-12,
            abs(public$weighted_right_fraction - fraction[3L]) < 1e-12,
            abs(public$exact_fraction - mean(exact)) < 1e-12)
          spot <- data.frame(medium = medium, element = element, state = region,
            cell_size = size, n_total = n, passed = TRUE)
          spot_checks[[medium]] <- spot; spot_done <- TRUE
        }
      }
    }
    table <- rbindlist(rows)
    saveRDS(list(element = element, rows = table, spot = spot), checkpoint)
    if (i %% 10L == 0L || i == length(elements)) {
      cat(medium, i, "/", length(elements), "species audited\n"); flush.console()
    }
  }
  all <- rbindlist(lapply(seq_along(elements), function(i)
    readRDS(file.path(folder, "checkpoints", sprintf("%03d.rds", i)))$rows))
  stopifnot(nrow(all) == nrow(original))
  fwrite(all, file.path(folder, "weighted_censoring_support.csv"))
  cat("Completed", medium, nrow(all), "cell/element/state rows\n"); flush.console()
  rm(chemistry, df, metadata, original, all); invisible(gc())
}
stopifnot(identical(tools::md5sum(sources), signature$sources),
          identical(tools::md5sum(input_files), signature$inputs))
write.csv(do.call(rbind, spot_checks), file.path(out, "public_api_spot_checks.csv"), row.names = FALSE)
writeLines(c("Completed all 80 soil and 81 sediment species.",
  "All original state/cell counts and Kish weight sizes reconciled.",
  "Source/input fingerprints unchanged. No concentration estimates or intervals refitted."),
  file.path(out, "COMPLETED.txt"))
