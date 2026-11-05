#!/usr/bin/env Rscript
# Reproduce the GitHub sediment examples with separately downloaded NGDB exports.
# Requires geoCensor 1.2.0 and sf; no input data are included in this repository.

select_records <- function(path, required, select, chunk_size = 10000L) {
  plain <- function(x) tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
  types <- readr::cols(.default = readr::col_character())
  header <- plain(readr::read_csv(path, n_max = 0, col_types = types,
    na = c("", "NA"), trim_ws = FALSE, name_repair = "check_unique",
    show_col_types = FALSE, progress = FALSE))
  if (!all(required %in% names(header))) stop("Missing columns in ", basename(path), ": ",
    paste(setdiff(required, names(header)), collapse = ", "), call. = FALSE)
  pieces <- indices <- list(); n_input <- 0
  callback <- readr::SideEffectChunkCallback$new(function(x, pos) {
    if (!identical(names(x), names(header)) || !all(vapply(x, is.character, logical(1))) ||
        nrow(readr::problems(x)) || pos != n_input + 1 || nrow(x) > chunk_size) {
      stop("CSV parsing/schema problem in ", basename(path), " near record ", n_input + 1, call. = FALSE)
    }
    keep <- which(select(x))
    if (length(keep)) {
      pieces[[length(pieces) + 1L]] <<- plain(x[keep, , drop = FALSE])
      indices[[length(indices) + 1L]] <<- n_input + keep
    }
    n_input <<- n_input + nrow(x)
  })
  withCallingHandlers(readr::read_csv_chunked(path, callback, chunk_size = chunk_size,
    col_types = types, na = c("", "NA"), trim_ws = FALSE,
    show_col_types = FALSE, progress = FALSE),
    warning = function(w) stop("CSV parsing warning: ", conditionMessage(w), call. = FALSE))
  list(data = if (length(pieces)) dplyr::bind_rows(pieces) else header,
       records = as.numeric(unlist(indices, use.names = FALSE)), scanned = n_input)
}

showcase_theme <- function() {
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold", colour = "#142D43"),
      plot.subtitle = ggplot2::element_text(size = 11, colour = "#41576A", margin = ggplot2::margin(b = 15)),
      plot.caption = ggplot2::element_text(size = 9, colour = "#41576A", hjust = 0, margin = ggplot2::margin(t = 14)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E4EAF0"),
      plot.margin = ggplot2::margin(18, 24, 18, 18), legend.position = "bottom")
}

save_showcase_figures <- function(case, summary, sensitivity, bounds, directory) {
  save_plot <- function(name, plot, height = 5.5) ggplot2::ggsave(
    file.path(directory, paste0(name, ".png")), plot, width = 10, height = height, dpi = 160, bg = "white")
  title <- if (case == "illinois") "Upper Illinois sediment records" else "Wisconsin Lead-Zinc project"
  caveat <- if (case == "illinois") "Sampling fractions vary; this is a project description, not a statewide estimate." else
    "Semiquantitative method flags retained; reporting classes and source attribution are unresolved."
  lead <- sensitivity[sensitivity$element == "Pb" & is.finite(sensitivity$geomean) & sensitivity$geomean > 0, ]
  if (nrow(lead)) {
    lead$mg_kg <- lead$geomean / 1000
    lead$cell_km <- lead$cell_size / 1000
    baseline <- unique(lead$unweighted_geomean / 1000)
    stopifnot(length(baseline) == 1L)
    p <- ggplot2::ggplot(lead, ggplot2::aes(x = cell_km, y = mg_kg)) +
      ggplot2::geom_hline(yintercept = baseline, colour = "#8B6345", linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("text", x = 12.5, y = baseline,
        label = paste0("Unweighted: ", format(round(baseline, 1), nsmall = 1), " mg/kg"),
        vjust = -0.8, colour = "#8B6345", size = 3.6) +
      ggplot2::geom_line(colour = "#087F8C", linewidth = 1.2) +
      ggplot2::geom_point(shape = 21, size = 4, stroke = 1.2, fill = "white", colour = "#087F8C") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", mg_kg)), vjust = 2, colour = "#142D43", size = 4) +
      ggplot2::scale_x_continuous(breaks = c(5, 10, 20), limits = c(3, 22)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.28))) +
      ggplot2::labs(title = paste0(title, ": lead"),
        subtitle = paste0("Descriptive geometric means | ", unique(lead$n_total),
          " usable records | 100 grid placements at each cell width"),
        x = "Declustering cell width (km)", y = "Lead geometric mean (mg/kg)",
        caption = paste("Teal: declustered fit. Dashed: unweighted fit on the same usable records.",
          "No spatial confidence intervals. Declustering does not establish concentrations in unsampled areas.", caveat, sep = "\n")) + showcase_theme()
    save_plot("lead_declustering", p)
  }
  reporting <- dplyr::bind_rows(lapply(c("n_exact", "n_left", "n_right"), function(name) {
    data.frame(element = summary$element, n_total = summary$n_total,
      fraction = ifelse(summary$n_total > 0, summary[[name]] / summary$n_total, NA_real_),
      type = c(n_exact = "Uncensored code", n_left = "Below a limit", n_right = "Above a limit")[[name]])
  }))
  reporting$element <- factor(reporting$element, levels = rev(c("Pb", "As", "Cd", "Cu", "Hg", "Zn")))
  labels <- summary
  labels$element <- factor(labels$element, levels = levels(reporting$element))
  labels$label <- ifelse(labels$n_total > 0, paste0("n = ", labels$n_total), "No records")
  p <- ggplot2::ggplot(reporting, ggplot2::aes(x = element, y = fraction, fill = type)) +
    ggplot2::geom_col(width = 0.62, na.rm = TRUE) + ggplot2::coord_flip(clip = "off") +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = element, y = 1.025, label = label),
      inherit.aes = FALSE, hjust = 0, size = 3.4, colour = "#41576A") +
    ggplot2::scale_fill_manual(values = c("Uncensored code" = "#087F8C", "Below a limit" = "#E8B35B", "Above a limit" = "#A574AC"),
      breaks = c("Uncensored code", "Below a limit", "Above a limit"), drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) paste0(100*x, "%"), limits = c(0, 1.2)) +
    ggplot2::labs(title = paste0(title, ": what was reported"),
      subtitle = "Counts describe reporting codes; uncensored coding does not certify measurement precision.",
      x = NULL, y = "Share of usable chemistry records", fill = NULL,
      caption = paste("Nondetects retain their limits. Missing concentrations are not replaced with zero or a fraction of a limit.", caveat, sep = "\n")) + showcase_theme()
  save_plot("reporting", p)
  shown <- sensitivity[sensitivity$element %in% summary$element, ]
  shown$relative_geomean <- shown$geomean / shown$unweighted_geomean
  shown$cell_km <- shown$cell_size / 1000
  p <- ggplot2::ggplot(shown, ggplot2::aes(x = cell_km, y = relative_geomean)) +
    ggplot2::geom_hline(yintercept = 1, colour = "#8B6345", linetype = "dashed") +
    ggplot2::geom_line(colour = "#087F8C", linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(colour = "#087F8C", size = 2.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~element, ncol = 3) +
    ggplot2::scale_x_continuous(breaks = c(5, 10, 20)) +
    ggplot2::labs(title = paste0(title, ": spatial sensitivity"),
      subtitle = "Each element is compared with its own unweighted fit on the same usable records.",
      x = "Cell width (km)", y = "Weighted / unweighted geometric mean",
      caption = paste("A ratio of 1 means no change. Missing fits remain missing; blank panels do not mean zero concentration.",
        "No confidence limits or comparisons of environmental risk are shown.", caveat, sep = "\n")) + showcase_theme()
  save_plot("metals_sensitivity", p, height = 6.5)
  if (nrow(bounds)) {
    selected <- bounds[bounds$element %in% c("Pb", "Cd"), ]
    p <- ggplot2::ggplot(selected, ggplot2::aes(x = threshold / 1000)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = cdf_lower, ymax = cdf_upper), fill = "#E8B35B", alpha = 0.45) +
      ggplot2::geom_line(ggplot2::aes(y = cdf_lower), colour = "#087F8C", linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = cdf_upper), colour = "#087F8C", linetype = "dashed", linewidth = 0.8) +
      ggplot2::scale_x_log10() + ggplot2::facet_wrap(~element, nrow = 1) +
      ggplot2::labs(title = paste0(title, ": information retained by limits"),
        subtitle = "Shading shows the possible empirical cumulative fractions given the reported values and limits.",
        x = "Reported concentration threshold (mg/kg; logarithmic scale)", y = "Empirical cumulative fraction",
        caption = paste("Equal record weights; all usable records. These are reporting-information bounds, not confidence bands.",
          "Bounds condition on encoded limits; no semiquantitative reporting-class ranges are inferred.", caveat, sep = "\n")) + showcase_theme()
    save_plot("reporting_bounds", p)
  }
}

run_showcase <- function(args) {
  usage <- "Usage: Rscript docs/examples/sediment-metals.R main.csv bestvalue.csv NEW_output_directory"
  if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  if (length(args) != 3L || any(startsWith(args, "--"))) stop(usage, call. = FALSE)
  suppressPackageStartupMessages(library(geoCensor))
  helper <- system.file("examples", "workflow_requirements.R", package = "geoCensor")
  if (!nzchar(helper)) stop("Install the final geoCensor 1.2.0 release and restart R", call. = FALSE)
  sys.source(helper, envir = environment()); geocensor_require_workflow("illinois")
  if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf to prepare EPSG:5070 coordinates for this example", call. = FALSE)
  target <- args[3]
  check_target <- function() {
    link <- Sys.readlink(target)
    if (file.exists(target) || (!is.na(link) && nzchar(link))) stop("Output already exists; choose a new directory", call. = FALSE)
    if (!dir.exists(dirname(target))) stop("Output parent directory must exist", call. = FALSE)
  }
  check_target()
  files <- c(metadata = args[1], chemistry = args[2])
  if (any(!file.exists(files) | dir.exists(files))) stop("Inputs must be readable CSV files", call. = FALSE)
  before <- tools::md5sum(files)
  cases <- data.frame(case = c("illinois", "wisconsin"), state = c("IL", "WI"),
    project = c("NAWQA Upper Illinois R. Basin", "Wisconsin Lead-Zinc"), stringsAsFactors = FALSE)
  main <- select_records(files[["metadata"]],
    c("lab_id", "state", "project_name", "primary_class", "long_wgs84", "lat_wgs84"), function(x) {
      key <- paste(toupper(trimws(x$state)), trimws(x$project_name))
      key %in% paste(cases$state, cases$project)
    })
  ids <- trimws(main$data$lab_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Selected metadata IDs must be nonblank and unique", call. = FALSE)
  if (any(c("x_m", "y_m") %in% names(main$data))) stop("Supply native geographic metadata, not previously projected columns", call. = FALSE)
  chemistry <- select_records(files[["chemistry"]], c("lab_id", "species", "unit", "qvalue"),
    function(x) trimws(x$lab_id) %in% ids)
  results <- list(); elements <- c("Pb", "As", "Cd", "Cu", "Hg", "Zn")
  for (i in seq_len(nrow(cases))) {
    config <- cases[i, ]
    keep <- which(toupper(trimws(main$data$state)) %in% config$state & trimws(main$data$project_name) %in% config$project)
    metadata_original <- main$data[keep, ]; metadata <- metadata_original
    if (!nrow(metadata)) stop("No records for ", config$project, call. = FALSE)
    if (anyNA(metadata$primary_class) || any(tolower(trimws(metadata$primary_class)) != "sediment")) {
      stop("Unexpected sampling medium; review the cohort before running the sediment example", call. = FALSE)
    }
    rows <- which(trimws(chemistry$data$lab_id) %in% trimws(metadata$lab_id))
    chem <- chemistry$data[rows, ]
    if (!nrow(chem)) stop("No chemistry for ", config$project, call. = FALSE)
    number <- function(x) {
      value <- readr::parse_double(trimws(x), na = c("", "NA"))
      if (nrow(readr::problems(value)) || any(!is.na(value) & !is.finite(value))) stop("Invalid coordinate text", call. = FALSE)
      value
    }
    lon <- number(metadata$long_wgs84); lat <- number(metadata$lat_wgs84)
    if (any(abs(lon) > 180 | abs(lat) > 90, na.rm = TRUE)) stop("Coordinates outside geographic bounds", call. = FALSE)
    located <- is.finite(lon) & is.finite(lat)
    metadata$x_m <- metadata$y_m <- rep(NA_real_, nrow(metadata))
    if (any(located)) {
      points <- sf::st_as_sf(data.frame(lon = lon[located], lat = lat[located]), coords = c("lon", "lat"), crs = 4326)
      xy <- sf::st_coordinates(sf::st_transform(points, 5070))
      metadata$x_m[located] <- xy[, 1]; metadata$y_m[located] <- xy[, 2]
    }
    result <- analyze_elements(chem, main_df = metadata, group_cols = "state",
      x_col = "x_m", y_col = "y_m", cell_sizes = c(5000, 10000, 20000), n_offsets = 10,
      stratify_methods = TRUE, audit_method_overlap = TRUE, audit_source_precision = TRUE)
    support <- sampling_support_audit(metadata, group_cols = "state",
      support_cols = intersect(c("depth", "sample_source", "method_collected", "primary_class",
        "specific_name", "addl_attr", "prep", "mesh_size", "collection_protocol"), names(metadata)),
      x_col = "x_m", y_col = "y_m", include_memberships = TRUE)
    summary <- dplyr::left_join(data.frame(element = elements), result$estimates, by = "element")
    absent <- is.na(summary$n_input)
    for (name in c("n_input", "n_total", "n_exact", "n_left", "n_right")) summary[[name]][absent] <- 0L
    summary$reporting_status <- ifelse(absent, "no_records", "reported")
    if (any(!is.na(summary$value_unit) & summary$value_unit != "ug/kg")) stop("Unexpected estimator units", call. = FALSE)
    summary$geomean_mg_kg <- summary$geomean / 1000
    bounds <- dplyr::bind_rows(lapply(elements[!absent], function(element) {
      value <- censored_cdf_bounds(build_element(chem, element), thresholds = 10^seq(0, 9, length.out = 181))$cdf
      cbind(element = element, value)
    }))
    if (any(!is.na(result$declustering$ci_low)) || any(!is.na(result$declustering$ci_high))) stop("Unexpected routine weighted confidence limits", call. = FALSE)
    counts <- data.frame(item = c("metadata_records", "chemistry_records", "species", "records_with_geographic_coordinates",
      "recorded_geographic_positions", "collection_dates_present"), value = c(nrow(metadata), nrow(chem), length(unique(chem$species)),
      sum(located), nrow(unique(data.frame(lon = lon[located], lat = lat[located]))),
      if ("date_collect" %in% names(metadata)) sum(!is.na(metadata$date_collect) & nzchar(trimws(metadata$date_collect))) else NA_integer_))
    results[[config$case]] <- list(metadata_original = metadata_original, metadata = metadata, chemistry = chem,
      metadata_index = main$records[keep], chemistry_index = chemistry$records[rows],
      result = result, support = support, summary = summary, bounds = bounds, inventory = counts)
  }
  if (!identical(before, tools::md5sum(files))) stop("Input files changed; no output written", call. = FALSE)
  check_target()
  if (!dir.create(target, showWarnings = FALSE)) stop("Could not create output", call. = FALSE)
  for (case in names(results)) {
    data <- results[[case]]; directory <- file.path(target, case); dir.create(directory)
    saveRDS(data, file.path(directory, "analysis.rds"))
    for (name in c("metadata_original", "metadata", "chemistry", "summary", "bounds", "inventory")) {
      readr::write_csv(data[[name]], file.path(directory, paste0(name, ".csv")))
    }
    readr::write_csv(data.frame(source_data_record = data$metadata_index), file.path(directory, "metadata_record_index.csv"))
    readr::write_csv(data.frame(source_data_record = data$chemistry_index), file.path(directory, "chemistry_record_index.csv"))
    for (name in names(data$result)) if (ncol(data$result[[name]])) readr::write_csv(data$result[[name]], file.path(directory, paste0(name, ".csv")))
    for (name in names(data$support)) if (is.data.frame(data$support[[name]]) && ncol(data$support[[name]])) {
      readr::write_csv(data$support[[name]], file.path(directory, paste0("support_", name, ".csv")))
    }
    save_showcase_figures(case, data$summary, data$result$declustering, data$bounds, directory)
  }
  provenance <- list(package_version = as.character(utils::packageVersion("geoCensor")),
    selection = cases, input_paths = normalizePath(files), input_md5 = before,
    scanned_metadata_records = main$scanned, scanned_chemistry_records = chemistry$scanned,
    projection = "Published EPSG:4326 WGS84 fields transformed to EPSG:5070; original locations not corrected",
    sf_version = as.character(utils::packageVersion("sf")), projection_versions = sf::sf_extSoftVersion(),
    cell_sizes_m = c(5000, 10000, 20000), grid_placements = 100,
    statistical_status = "Descriptive; no routine weighted confidence limits; no source attribution or risk classification")
  saveRDS(provenance, file.path(target, "provenance.rds"))
  writeLines(capture.output(sessionInfo()), file.path(target, "session_info.txt"))
  writeLines(c("Separate metadata-selected sediment project cohorts; no selection on concentrations.",
    "All reported species retained; six elements displayed. Missing observations and metadata remain missing.",
    "Uncensored coding does not establish analytical precision. No semiquantitative class bounds are inferred.",
    "No statewide background, contamination attribution, exposure history, health-risk or ore classification.",
    "Use a new output directory for every run. Partial write-failure outputs are retained."), file.path(target, "README.txt"))
  paths <- list.files(target, recursive = TRUE, full.names = TRUE)
  readr::write_csv(data.frame(file = substring(paths, nchar(target) + 2L), bytes = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths))), file.path(target, "output_checksums.csv"))
  message("Separate sediment examples saved to ", normalizePath(target))
  invisible(results)
}

if (sys.nframe() == 0L) run_showcase(commandArgs(trailingOnly = TRUE))
