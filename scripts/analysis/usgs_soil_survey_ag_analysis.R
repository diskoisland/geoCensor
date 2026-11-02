#!/usr/bin/env Rscript

# USGS Data Series 801, Appendix 3A: national A-horizon soil silver summary.
#
# Workbook structure (Appendix_3a_Ahorizon_18Sept2013.xlsx):
#   * The script finds the analytical worksheet and header row by searching for
#     an Ag field (`A_Ag` in the supplied workbook), rather than hard-coding the
#     worksheet name or row number.
#   * The supplied workbook has one analytical worksheet (`A_Horizon`), field
#     names on row 13, and units on the following row. `A_Ag` is column 36 and
#     is reported in mg/kg.
#
# Censoring and missingness:
#   * Unlike the heterogeneous NGDB files, Appendix 3A stores the qualifier and
#     reporting limit together in the Ag cell (for example, `<1`). There is no
#     separate Ag qualifier, reporting-limit, MDL, or LOD column.
#   * `<x`/`<=x`/`≤x` is mapped to the interval [0, x], and `>x`/`>=x`/`≥x`
#     is mapped to [x, Inf]. Numeric positive values are exact observations.
#   * `N.S.` (no sample) and `INS` (insufficient sample) are excluded because
#     they are not analytical non-detects. A bare `N.D.` cannot be assigned an
#     interval without a reporting limit and therefore triggers an error.
#
# Statistical method:
#   The interval-censored lognormal MLE is not reimplemented here. The script
#   converts Appendix 3A into geoCensor's standard censoring columns and calls
#   the existing `fit_cen_lognorm()` helper from R/censoring.R, which uses
#   survival::Surv(..., type = "interval2") and survival::survreg().
#
# Arithmetic mean, median, geometric mean, minimum, and maximum describe the
# detected concentrations only; they are labelled accordingly in the output.
# With 99%+ censoring at one reporting limit, these descriptive statistics and
# the fitted MLE answer different questions and should not be interchanged.

suppressPackageStartupMessages({
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required.")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required by fit_cen_lognorm().")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required by fit_cen_lognorm().")
  }
})

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_start <- if (length(file_arg)) {
  # Some Rscript launchers encode spaces in --file paths as "~+~".
  script_file <- gsub("~[+]~", " ", sub("^--file=", "", file_arg[1]))
  dirname(normalizePath(script_file, mustWork = TRUE))
} else {
  getwd()
}

find_repo_root <- function(start = script_start) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "R", "censoring.R"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root containing DESCRIPTION and R/censoring.R.")
    }
    path <- parent
  }
}

find_appendix_workbook <- function(repo_root) {
  search_dirs <- unique(c(
    file.path(repo_root, "data", "master", "USGS_soil_Survey"),
    file.path(repo_root, "data", "master", "USGS_soil_survey"),
    file.path(repo_root, "..", "..", "USGS_soil_Survey"),
    file.path(repo_root, "..", "..", "USGS_soil_survey")
  ))
  search_dirs <- search_dirs[dir.exists(search_dirs)]
  files <- unlist(lapply(
    search_dirs,
    list.files,
    pattern = "^Appendix_3a.*[.]xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  ), use.names = FALSE)
  files <- files[!grepl("^~[$]", basename(files))]
  files <- unique(normalizePath(files, mustWork = TRUE))
  if (length(files) != 1L) {
    stop("Expected exactly one Appendix 3A workbook; found ", length(files),
         if (length(files)) paste0(": ", paste(files, collapse = ", ")) else ".")
  }
  files
}

find_ag_sheet <- function(workbook, max_header_rows = 100L) {
  sheets <- readxl::excel_sheets(workbook)
  hits <- lapply(sheets, function(sheet) {
    preview <- suppressMessages(readxl::read_excel(
      workbook, sheet = sheet, col_names = FALSE, col_types = "text",
      n_max = max_header_rows
    ))
    values <- as.matrix(preview)
    is_ag <- matrix(
      trimws(as.character(values)) %in% c("A_Ag", "Ag"),
      nrow = nrow(values), ncol = ncol(values)
    )
    locations <- which(is_ag, arr.ind = TRUE)
    if (!nrow(locations)) return(NULL)
    data.frame(
      sheet = sheet,
      header_row = locations[, "row"],
      ag_column = locations[, "col"],
      ag_name = values[locations],
      stringsAsFactors = FALSE
    )
  })
  hits <- do.call(rbind, hits[!vapply(hits, is.null, logical(1))])
  if (is.null(hits) || nrow(hits) != 1L) {
    stop("Expected one Ag analytical column across all worksheets; found ",
         if (is.null(hits)) 0L else nrow(hits), ".")
  }
  hits
}

parse_ag_cells <- function(values) {
  raw <- trimws(as.character(values))
  raw[raw == ""] <- NA_character_
  numeric_value <- suppressWarnings(as.numeric(raw))

  left <- grepl("^(<|<=|≤)[[:space:]]*[0-9]", raw)
  right <- grepl("^(>|>=|≥)[[:space:]]*[0-9]", raw)
  strip_qualifier <- function(x) {
    suppressWarnings(as.numeric(sub("^(<=|>=|<|>|≤|≥)[[:space:]]*", "", x)))
  }
  left_limit <- rep(NA_real_, length(raw))
  right_limit <- rep(NA_real_, length(raw))
  left_limit[left] <- strip_qualifier(raw[left])
  right_limit[right] <- strip_qualifier(raw[right])

  exact <- is.finite(numeric_value) & numeric_value > 0 & !left & !right
  missing_code <- toupper(raw) %in% c("N.S.", "NS", "INS") | is.na(raw)
  no_limit_nd <- toupper(raw) %in% c("N.D.", "ND")
  if (any(no_limit_nd, na.rm = TRUE)) {
    stop("Ag contains bare N.D. values with no reporting limit; an interval cannot be inferred.")
  }
  recognized <- exact | left | right | missing_code
  if (any(!recognized)) {
    stop("Unrecognized Ag values: ", paste(unique(raw[!recognized]), collapse = ", "))
  }
  if (any(left & (!is.finite(left_limit) | left_limit <= 0)) ||
      any(right & (!is.finite(right_limit) | right_limit <= 0))) {
    stop("All censored Ag observations must contain a positive numeric limit.")
  }

  data.frame(
    detected_val = ifelse(exact, numeric_value, NA_real_),
    left_cen = left,
    left_limit = left_limit,
    right_cen = right,
    right_limit = right_limit,
    stringsAsFactors = FALSE
  )
}

repo_root <- find_repo_root()
workbook <- find_appendix_workbook(repo_root)
location <- find_ag_sheet(workbook)

raw_data <- suppressMessages(readxl::read_excel(
  workbook,
  sheet = location$sheet,
  skip = location$header_row - 1L,
  col_types = "text"
))
if (!location$ag_name %in% names(raw_data)) {
  stop("Ag column was not preserved when the analytical table was imported.")
}

# The first post-header row contains units, not a sample observation.
unit <- trimws(raw_data[[location$ag_name]][1])
if (!identical(tolower(unit), "mg/kg")) {
  stop("Expected Ag units of mg/kg; found '", unit, "'.")
}
sample_data <- raw_data[-1, , drop = FALSE]
ag_intervals <- parse_ag_cells(sample_data[[location$ag_name]])
usable <- ag_intervals$left_cen | ag_intervals$right_cen |
  is.finite(ag_intervals$detected_val)
ag_intervals <- ag_intervals[usable, , drop = FALSE]

# Reuse the project's published NGDB censoring-aware implementation.
source(file.path(repo_root, "R", "distribution_diagnostics.R"), local = TRUE)
source(file.path(repo_root, "R", "flag_analytical_methods.R"), local = TRUE)
source(file.path(repo_root, "R", "censoring.R"), local = TRUE)
mle <- fit_cen_lognorm(ag_intervals)

detections <- ag_intervals$detected_val[is.finite(ag_intervals$detected_val)]
n_samples <- nrow(ag_intervals)
n_detections <- length(detections)
if (!n_detections) stop("No detected Ag concentrations were found.")

summary_table <- data.frame(
  dataset = "USGS Data Series 801 Appendix 3A A-horizon soils",
  worksheet = location$sheet,
  ag_column = location$ag_name,
  units = "mg/kg",
  n_samples = n_samples,
  n_detections = n_detections,
  detection_frequency = n_detections / n_samples,
  arithmetic_mean_detections = mean(detections),
  median_detections = stats::median(detections),
  geometric_mean_detections = exp(mean(log(detections))),
  mle_geometric_mean = mle$geomean,
  mle_ci_95_lower = mle$ci_low,
  mle_ci_95_upper = mle$ci_high,
  minimum_detected = min(detections),
  maximum_detected = max(detections),
  mle_method = mle$method,
  stringsAsFactors = FALSE
)

results_dir <- file.path(repo_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(results_dir, "usgs_soil_survey_ag_summary.csv")
utils::write.csv(summary_table, output_file, row.names = FALSE, na = "")

cat("USGS Soil Survey (Appendix 3A)\n\n")
cat(sprintf("Samples: %d\n", n_samples))
cat(sprintf("Detections: %d\n", n_detections))
cat(sprintf("Detection frequency: %.2f%%\n", 100 * n_detections / n_samples))
cat(sprintf("Arithmetic mean (detections): %.6g mg/kg\n", mean(detections)))
cat(sprintf("Median (detections): %.6g mg/kg\n", stats::median(detections)))
cat(sprintf("Geometric mean (detections): %.6g mg/kg\n", exp(mean(log(detections)))))
cat(sprintf("MLE geometric mean: %.6g mg/kg\n", mle$geomean))
cat(sprintf("95%% CI: %.6g to %.6g mg/kg\n", mle$ci_low, mle$ci_high))
cat(sprintf("Minimum detected Ag: %.6g mg/kg\n", min(detections)))
cat(sprintf("Maximum detected Ag: %.6g mg/kg\n", max(detections)))
cat(sprintf("\nSummary written to: %s\n", output_file))
