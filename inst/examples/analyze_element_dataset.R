#!/usr/bin/env Rscript

# Usage: Rscript analyze_element_dataset.R bestvalue.csv Ag output/Ag [--stratify-methods] [--audit-source-precision]
# Uses the installed geoCensor package. All concentrations are in ug/kg.
args <- commandArgs(trailingOnly = TRUE)
stratify_methods <- any(args == "--stratify-methods")
if (sum(args == "--stratify-methods") > 1L) stop("Supply --stratify-methods only once", call. = FALSE)
args <- args[args != "--stratify-methods"]
audit_source_precision <- any(args == "--audit-source-precision")
if (sum(args == "--audit-source-precision") > 1L) stop("Supply --audit-source-precision only once", call. = FALSE)
args <- args[args != "--audit-source-precision"]
if (length(args) != 3L) {
  stop("Usage: Rscript analyze_element_dataset.R chemistry.csv element output_dir [--stratify-methods] [--audit-source-precision]",
       call. = FALSE)
}
suppressPackageStartupMessages(library(geoCensor))
required_version <- "1.2.0.9006"
if (utils::packageVersion("geoCensor") < required_version) {
  stop("Install geoCensor ", required_version, " or newer to use the requested workflow.", call. = FALSE)
}
output_dir <- args[3]
check_output_directory <- function() {
  if (file.exists(output_dir) && !dir.exists(output_dir)) stop("Output path is not a directory", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; choose a new directory to preserve previous results", call. = FALSE)
  }
}
check_output_directory()

chemistry <- read_usgs_chemistry(args[1])
required <- c("lab_id", "species", "qvalue", "unit")
if (!all(required %in% names(chemistry))) {
  stop("Missing chemistry columns: ",
       paste(setdiff(required, names(chemistry)), collapse = ", "))
}
observations <- build_element(chemistry, element = args[2], to_ugkg = TRUE)
if (!nrow(observations)) stop("No observations for element: ", args[2])
if (any(!is.finite(observations$unit_factor))) {
  stop("Unsupported or missing concentration units; convert to supported mass-fraction units first")
}
ids <- trimws(as.character(observations$lab_id))
if (anyNA(ids) || any(ids == "")) stop("Missing or blank lab_id values")
if (anyDuplicated(ids)) stop("Duplicate lab_id values; select one result per sample first")
audit <- dplyr::count(observations, .data$validation_status, name = "n_observations")
batch_args <- list(chem_df = chemistry, elements = args[2], stratify_methods = stratify_methods)
if (audit_source_precision) batch_args$audit_source_precision <- TRUE
batch <- do.call(analyze_elements, batch_args)
estimate <- batch$estimates
comparison <- batch$distributions
methods <- batch$methods
plot <- if (is.finite(estimate$geomean)) {
  plot_lognormal_fit(
    observations, title = paste(args[2], "censored lognormal fit"),
    show_nonparametric = TRUE
  )
} else NULL

check_output_directory()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory", call. = FALSE)
readr::write_csv(observations, file.path(output_dir, "observations.csv"))
readr::write_csv(audit, file.path(output_dir, "validation_summary.csv"))
readr::write_csv(methods, file.path(output_dir, "methods.csv"))
readr::write_csv(estimate, file.path(output_dir, "pooled_estimate.csv"))
readr::write_csv(comparison, file.path(output_dir, "distribution_comparison.csv"))
if (stratify_methods) readr::write_csv(batch$method_estimates, file.path(output_dir, "method_estimates.csv"))
if (audit_source_precision) {
  readr::write_csv(batch$source_precision, file.path(output_dir, "source_precision.csv"))
  readr::write_csv(batch$source_precision_summary, file.path(output_dir, "source_precision_summary.csv"))
}
if (!is.null(plot)) {
  ggplot2::ggsave(file.path(output_dir, "lognormal_fit.png"), plot,
                 width = 8, height = 5, dpi = 180)
}
analysis <- list(observations = observations, validation = audit, estimate = estimate,
                 comparison = comparison, methods = methods)
if (stratify_methods) analysis$method_estimates <- batch$method_estimates
if (audit_source_precision) {
  analysis$source_precision <- batch$source_precision
  analysis$source_precision_summary <- batch$source_precision_summary
}
saveRDS(analysis, file.path(output_dir, "analysis.rds"))
writeLines(c(
  paste("Input:", normalizePath(args[1])),
  "Chemistry import: text fields; blank/literal NA missing; malformed CSV rejected.",
  paste("Element:", args[2]),
  "Concentration unit: ug/kg",
  paste("Input rows for element:", nrow(observations)),
  paste("Usable observations:", estimate$n_total),
  paste("Excluded observations:", nrow(observations) - estimate$n_total),
  paste("Estimation method:", estimate$method),
  paste("Method review:", estimate$method_review_status),
  paste("Unweighted method stratification:", stratify_methods),
  paste("Source-summary precision audit:", audit_source_precision),
  if (audit_source_precision) "Source-precision candidates/deltas retain input units and do not replace estimator inputs.",
  "Method review does not infer reporting-class bounds or validate measurement precision."
), file.path(output_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Analysis saved to ", normalizePath(output_dir))
