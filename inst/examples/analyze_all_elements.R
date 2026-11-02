#!/usr/bin/env Rscript
# Analyze all species; run separately for soil and sediment datasets.
# Optional metadata can supply projected coordinates and analysis groups.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L, 6L, 7L)) {
  stop(paste(
    "Usage: Rscript analyze_all_elements.R chemistry.csv output_dir",
    "[metadata.csv [x_col y_col cell_sizes [group_cols]]]",
    "Separate cell sizes and grouping columns with commas.",
    "Example: bestvalue.csv output/all main_projected.csv x_m y_m 5000,10000,20000 state",
    sep = "\n"
  ), call. = FALSE)
}
suppressPackageStartupMessages(library(geoCensor))
read_input <- function(path) readr::read_csv(
  path, col_types = readr::cols(lab_id = readr::col_character()), show_col_types = FALSE
)
chemistry <- read_input(args[1])
metadata <- if (length(args) >= 3L) read_input(args[3]) else NULL
sizes <- if (length(args) >= 6L) as.numeric(strsplit(args[6], ",", fixed = TRUE)[[1]]) else NULL
groups <- if (length(args) == 7L) strsplit(args[7], ",", fixed = TRUE)[[1]] else character()
result <- analyze_elements(
  chemistry, main_df = metadata, group_cols = groups, cell_sizes = sizes,
  x_col = if (length(args) >= 6L) args[4] else NULL,
  y_col = if (length(args) >= 6L) args[5] else NULL
)
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (name in names(result)) {
  if (ncol(result[[name]])) readr::write_csv(result[[name]], file.path(output_dir, paste0(name, ".csv")))
}
saveRDS(result, file.path(output_dir, "analysis.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
writeLines(c(
  paste("Input:", normalizePath(args[1])),
  paste("Elements:", length(unique(result$estimates$element))),
  paste("Analysis groups:", if (length(groups)) paste(groups, collapse = ", ") else "pooled"),
  paste("Declustering:", if (is.null(sizes)) "not requested" else paste(sizes, collapse = ", ")),
  "Concentration unit: ug/kg",
  "Weighted confidence limits are not computed."
), file.path(output_dir, "run_summary.txt"))
message("All-element analysis saved to ", normalizePath(output_dir))
