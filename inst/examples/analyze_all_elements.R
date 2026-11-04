#!/usr/bin/env Rscript
# Analyze all species; run separately for soil and sediment datasets.
# Optional metadata can supply projected coordinates and analysis groups.
args <- commandArgs(trailingOnly = TRUE)
stratify_methods <- any(args == "--stratify-methods")
if (sum(args == "--stratify-methods") > 1L) stop("Supply --stratify-methods only once", call. = FALSE)
args <- args[args != "--stratify-methods"]
audit_method_overlap <- any(args == "--method-overlap")
if (sum(args == "--method-overlap") > 1L) stop("Supply --method-overlap only once", call. = FALSE)
args <- args[args != "--method-overlap"]
audit_source_precision <- any(args == "--audit-source-precision")
if (sum(args == "--audit-source-precision") > 1L) stop("Supply --audit-source-precision only once", call. = FALSE)
args <- args[args != "--audit-source-precision"]
if (!length(args) %in% c(2L, 3L, 6L, 7L)) {
  stop(paste(
    "Usage: Rscript analyze_all_elements.R chemistry.csv output_dir",
    "[metadata.csv [x_col y_col cell_sizes [group_cols]]] [--stratify-methods] [--method-overlap] [--audit-source-precision]",
    "Separate cell sizes and grouping columns with commas.",
    "Example: bestvalue.csv output/all main_projected.csv x_m y_m 5000,10000,20000 state",
    sep = "\n"
  ), call. = FALSE)
}
if (audit_method_overlap && length(args) < 6L) {
  stop("--method-overlap requires metadata, projected x_col/y_col and cell_sizes", call. = FALSE)
}
suppressPackageStartupMessages(library(geoCensor))
required_version <- "1.2.0.9006"
if (utils::packageVersion("geoCensor") < required_version) {
  stop("Install geoCensor ", required_version, " or newer to use the requested workflow.", call. = FALSE)
}
output_dir <- args[2]
check_output_directory <- function() {
  if (file.exists(output_dir) && !dir.exists(output_dir)) stop("Output path is not a directory", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a new directory to preserve previous results", call. = FALSE)
  }
}
check_output_directory()
chemistry <- read_usgs_chemistry(args[1])
metadata <- NULL
if (length(args) >= 3L) {
  metadata <- readr::read_csv(args[3], col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"), trim_ws = FALSE, name_repair = "check_unique", show_col_types = FALSE, progress = FALSE)
  if (nrow(readr::problems(metadata))) stop("Metadata CSV parsing problems; resolve malformed records before analysis", call. = FALSE)
  if (length(args) >= 6L) {
    for (column in unique(args[4:5])) {
      if (!column %in% names(metadata)) stop("Missing projected coordinate column: ", column, call. = FALSE)
      text <- trimws(metadata[[column]])
      nonfinite <- !is.na(text) & toupper(text) %in% c("INF", "+INF", "-INF", "NAN")
      regular <- which(!nonfinite)
      parsed <- readr::parse_double(text[regular], na = c("", "NA"))
      problems <- readr::problems(parsed)
      if (nrow(problems)) stop("Invalid numeric coordinate column '", column, "' at input rows: ",
        paste(utils::head(regular[problems$row], 10L), collapse = ", "), call. = FALSE)
      value <- rep(NA_real_, length(text))
      value[regular] <- parsed
      lookup <- c(INF = Inf, "+INF" = Inf, "-INF" = -Inf, NAN = NaN)
      value[nonfinite] <- unname(lookup[toupper(text[nonfinite])])
      metadata[[column]] <- value
    }
  }
}
sizes <- if (length(args) >= 6L) as.numeric(strsplit(args[6], ",", fixed = TRUE)[[1]]) else NULL
groups <- if (length(args) == 7L) strsplit(args[7], ",", fixed = TRUE)[[1]] else character()
batch_args <- list(
  chem_df = chemistry, main_df = metadata, group_cols = groups, cell_sizes = sizes,
  x_col = if (length(args) >= 6L) args[4] else NULL,
  y_col = if (length(args) >= 6L) args[5] else NULL,
  stratify_methods = stratify_methods,
  audit_method_overlap = audit_method_overlap
)
if (audit_source_precision) batch_args$audit_source_precision <- TRUE
result <- do.call(analyze_elements, batch_args)
check_output_directory()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory", call. = FALSE)
for (name in names(result)) {
  if (ncol(result[[name]])) readr::write_csv(result[[name]], file.path(output_dir, paste0(name, ".csv")))
}
saveRDS(result, file.path(output_dir, "analysis.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
writeLines(c(
  paste("Input:", normalizePath(args[1])),
  "Chemistry import: text fields; blank/literal NA missing; malformed CSV rejected.",
  if (!is.null(metadata)) "Metadata import: text fields; only explicitly supplied projected coordinates parsed numeric.",
  paste("Elements:", length(unique(result$estimates$element))),
  paste("Analysis groups:", if (length(groups)) paste(groups, collapse = ", ") else "pooled"),
  paste("Declustering:", if (is.null(sizes)) "not requested" else paste(sizes, collapse = ", ")),
  "Concentration unit: ug/kg",
  "Weighted confidence limits are not computed.",
  paste("Unweighted method stratification:", stratify_methods),
  paste("Geographic method overlap:", audit_method_overlap),
  paste("Source-summary precision audit:", audit_source_precision),
  if (audit_source_precision) "Source-precision candidates/deltas retain input units and do not replace estimator inputs.",
  if (audit_method_overlap) "Method overlap uses common grids; placement ranges are not confidence limits or paired-sample evidence.",
  "methods.csv records input and unweighted fit-subset method counts; unresolved flags stay NA."
), file.path(output_dir, "run_summary.txt"))
message("All-element analysis saved to ", normalizePath(output_dir))
