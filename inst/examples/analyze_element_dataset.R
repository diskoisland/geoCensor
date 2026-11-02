#!/usr/bin/env Rscript

# Usage: Rscript analyze_element_dataset.R bestvalue.csv Ag output/Ag
# Uses the installed geoCensor package. All concentrations are in ug/kg.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript analyze_element_dataset.R chemistry.csv element output_dir",
       call. = FALSE)
}
suppressPackageStartupMessages(library(geoCensor))

chemistry <- readr::read_csv(args[1], col_types = readr::cols(lab_id = readr::col_character()),
                             show_col_types = FALSE)
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
estimate <- fit_cen_lognorm(observations)
comparison <- compare_censored_distributions(observations)
plot <- if (is.finite(estimate$geomean)) {
  plot_lognormal_fit(
    observations, title = paste(args[2], "censored lognormal fit"),
    show_nonparametric = TRUE
  )
} else NULL

output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(observations, file.path(output_dir, "observations.csv"))
readr::write_csv(audit, file.path(output_dir, "validation_summary.csv"))
readr::write_csv(estimate, file.path(output_dir, "pooled_estimate.csv"))
readr::write_csv(comparison, file.path(output_dir, "distribution_comparison.csv"))
if (!is.null(plot)) {
  ggplot2::ggsave(file.path(output_dir, "lognormal_fit.png"), plot,
                 width = 8, height = 5, dpi = 180)
}
saveRDS(list(observations = observations, validation = audit, estimate = estimate,
             comparison = comparison), file.path(output_dir, "analysis.rds"))
writeLines(c(
  paste("Input:", normalizePath(args[1])),
  paste("Element:", args[2]),
  "Concentration unit: ug/kg",
  paste("Input rows for element:", nrow(observations)),
  paste("Usable observations:", estimate$n_total),
  paste("Excluded observations:", nrow(observations) - estimate$n_total),
  paste("Estimation method:", estimate$method)
), file.path(output_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Analysis saved to ", normalizePath(output_dir))
