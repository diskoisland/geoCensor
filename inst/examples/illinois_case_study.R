#!/usr/bin/env Rscript
# A real-data workflow; input CSVs are downloaded separately from USGS.
run_illinois <- function(args) {
  usage <- paste(
    "Usage: Rscript illinois_case_study.R main.csv bestvalue.csv new_output_dir",
    "[--projected-metadata=metadata_with_x_m_y_m.csv] [--chunk-size=10000]",
    "Selects IL / NAWQA Upper Illinois R. Basin, without filtering chemistry.",
    "Coordinates: published long_wgs84/lat_wgs84, transformed to EPSG:5070.",
    "Alternatively supply matching sample IDs and x_m/y_m already in EPSG:5070.",
    "The sf package is needed only when preparing projected coordinates here.",
    sep = "\n")
  if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  projected <- NULL; chunk_size <- 10000L; positional <- character(); seen <- character()
  for (arg in args) {
    if (!startsWith(arg, "--")) { positional <- c(positional, arg); next }
    token <- regmatches(arg, regexec("^--(projected-metadata|chunk-size)=(.+)$", arg))[[1]]
    if (length(token) != 3L || token[2] %in% seen) stop(usage, call. = FALSE)
    seen <- c(seen, token[2])
    if (token[2] == "projected-metadata") projected <- token[3] else {
      if (!grepl("^[0-9]+$", token[3])) stop("chunk-size must be a positive integer", call. = FALSE)
      chunk_size <- as.numeric(token[3])
      if (!is.finite(chunk_size) || chunk_size < 1 || chunk_size > .Machine$integer.max) {
        stop("chunk-size must be a positive integer", call. = FALSE)
      }
    }
  }
  if (length(positional) != 3L) stop(usage, call. = FALSE)
  suppressPackageStartupMessages(library(geoCensor))
  helper <- system.file("examples", "workflow_requirements.R", package = "geoCensor")
  if (!nzchar(helper)) stop("Install the final geoCensor 1.2.0 release and restart R", call. = FALSE)
  sys.source(helper, envir = environment())
  geocensor_require_workflow("illinois")
  target <- positional[3]
  check_target <- function() {
    link <- Sys.readlink(target)
    if (file.exists(target) || (!is.na(link) && nzchar(link))) stop("Output already exists; choose a new directory", call. = FALSE)
    if (!dir.exists(dirname(target))) stop("Output parent directory must exist", call. = FALSE)
  }
  check_target()
  files <- c(metadata = positional[1], chemistry = positional[2],
    if (!is.null(projected)) c(projected_metadata = projected))
  if (any(!file.exists(files) | dir.exists(files))) stop("Inputs must be readable CSV files", call. = FALSE)
  before <- tools::md5sum(files)

  # Scan every batch, including unselected records, while retaining only the project.
  select_csv <- function(path, required, select) {
    types <- readr::cols(.default = readr::col_character())
    plain <- function(x) tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
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
      source_rows = as.numeric(unlist(indices, use.names = FALSE)), n_input = n_input)
  }
  project <- "NAWQA Upper Illinois R. Basin"
  main <- select_csv(files[["metadata"]], c("lab_id", "state", "project_name", "long_wgs84", "lat_wgs84"),
    function(x) toupper(trimws(x$state)) %in% "IL" & trimws(x$project_name) %in% project)
  metadata <- main$data; ids <- trimws(metadata$lab_id)
  if (!length(ids)) stop("No metadata rows match IL / ", project, call. = FALSE)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Project metadata IDs must be nonblank and unique", call. = FALSE)
  if (any(c("x_m", "y_m") %in% names(metadata))) stop("Use native main.csv; supply projected coordinates with --projected-metadata", call. = FALSE)
  chemistry <- select_csv(files[["chemistry"]], c("lab_id", "species", "unit", "qvalue"),
    function(x) trimws(x$lab_id) %in% ids)
  if (!nrow(chemistry$data)) stop("No chemistry matches the selected sample IDs", call. = FALSE)
  numeric_column <- function(x, name) {
    value <- readr::parse_double(trimws(x), na = c("", "NA"))
    if (nrow(readr::problems(value))) stop("Invalid numeric text in ", name, call. = FALSE)
    if (any(!is.na(value) & !is.finite(value))) stop("Nonfinite coordinates in ", name, call. = FALSE)
    value
  }
  lon <- numeric_column(metadata$long_wgs84, "long_wgs84")
  lat <- numeric_column(metadata$lat_wgs84, "lat_wgs84")
  if (any(abs(lon) > 180 | abs(lat) > 90, na.rm = TRUE)) stop("Geographic coordinates outside valid bounds", call. = FALSE)
  if (is.null(projected)) {
    if (!requireNamespace("sf", quietly = TRUE)) stop("Install sf for coordinate preparation, or supply --projected-metadata in EPSG:5070", call. = FALSE)
    located <- which(is.finite(lon) & is.finite(lat))
    metadata$x_m <- metadata$y_m <- rep(NA_real_, nrow(metadata))
    if (length(located)) {
      points <- sf::st_as_sf(data.frame(lon = lon[located], lat = lat[located]),
        coords = c("lon", "lat"), crs = 4326)
      xy <- sf::st_coordinates(sf::st_transform(points, 5070))
      metadata$x_m[located] <- xy[, 1]; metadata$y_m[located] <- xy[, 2]
    }
    projection <- list(method = "sf EPSG:4326 to EPSG:5070 on published WGS84 fields",
      sf_version = as.character(utils::packageVersion("sf")), external_versions = sf::sf_extSoftVersion())
  } else {
    supplied <- select_csv(projected, c("lab_id", "x_m", "y_m"), function(x) trimws(x$lab_id) %in% ids)
    keys <- trimws(supplied$data$lab_id)
    if (length(keys) != length(ids) || anyDuplicated(keys) || !setequal(keys, ids)) {
      stop("Projected metadata must supply exactly one row per selected sample ID", call. = FALSE)
    }
    position <- match(ids, keys)
    metadata$x_m <- numeric_column(supplied$data$x_m[position], "x_m")
    metadata$y_m <- numeric_column(supplied$data$y_m[position], "y_m")
    projection <- list(method = "User-supplied EPSG:5070 metre coordinates; CRS asserted by user",
      source_rows_in_sample_order = supplied$source_rows[position])
  }
  support_cols <- intersect(c("depth", "sample_source", "method_collected", "primary_class",
    "specific_name", "addl_attr", "prep", "mesh_size", "collection_protocol"), names(metadata))
  support <- sampling_support_audit(metadata, group_cols = "state", support_cols = support_cols,
    x_col = "x_m", y_col = "y_m", include_memberships = TRUE)
  result <- analyze_elements(chemistry$data, main_df = metadata, group_cols = "state",
    x_col = "x_m", y_col = "y_m", cell_sizes = c(5000, 10000, 20000), n_offsets = 10,
    stratify_methods = TRUE, audit_method_overlap = TRUE, audit_source_precision = TRUE)
  ag_bounds <- if ("Ag" %in% chemistry$data$species) censored_cdf_bounds(
    build_element(chemistry$data, "Ag"), thresholds = 10^seq(-1, 6, length.out = 121)) else NULL
  if (!identical(before, tools::md5sum(files))) stop("Input bytes changed during the run; no output written", call. = FALSE)
  check_target()
  if (!dir.create(target, recursive = FALSE, showWarnings = FALSE)) stop("Could not create a new output directory", call. = FALSE)
  saveRDS(result, file.path(target, "analysis.rds")); saveRDS(support, file.path(target, "support.rds"))
  if (!is.null(ag_bounds)) {
    saveRDS(ag_bounds, file.path(target, "Ag_cdf_bounds.rds"))
    readr::write_csv(ag_bounds$cdf, file.path(target, "Ag_cdf_bounds.csv"))
    ggplot2::ggsave(file.path(target, "Ag_cdf_bounds.png"), plot_censored_cdf_bounds(
      ag_bounds, title = "Upper Illinois Ag: empirical reporting bounds"), width = 8, height = 5, dpi = 160)
  }
  readr::write_csv(main$data, file.path(target, "selected_metadata_original.csv"))
  readr::write_csv(metadata, file.path(target, "metadata_projected.csv"))
  readr::write_csv(chemistry$data, file.path(target, "selected_chemistry.csv"))
  readr::write_csv(data.frame(lab_id = metadata$lab_id, source_data_record = main$source_rows),
    file.path(target, "metadata_record_index.csv"))
  readr::write_csv(data.frame(source_data_record = chemistry$source_rows), file.path(target, "chemistry_record_index.csv"))
  for (name in names(result)) if (ncol(result[[name]])) readr::write_csv(result[[name]], file.path(target, paste0(name, ".csv")))
  for (name in names(support)) if (is.data.frame(support[[name]]) && ncol(support[[name]])) {
    readr::write_csv(support[[name]], file.path(target, paste0("support_", name, ".csv")))
  }
  observed_positions <- unique(data.frame(longitude = lon, latitude = lat))
  located <- is.finite(observed_positions$longitude) & is.finite(observed_positions$latitude)
  counts <- data.frame(item = c("metadata_records", "chemistry_records", "species",
    "recorded_geographic_positions", "collection_dates_present"), value = c(nrow(metadata),
    nrow(chemistry$data), length(unique(chemistry$data$species)), sum(located),
    if ("date_collect" %in% names(metadata)) sum(!is.na(metadata$date_collect) & nzchar(trimws(metadata$date_collect))) else NA))
  readr::write_csv(counts, file.path(target, "inventory.csv"))
  map <- ggplot2::ggplot(observed_positions[located, ], ggplot2::aes(x = longitude, y = latitude)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.7) + ggplot2::coord_quickmap() +
    ggplot2::labs(title = "Upper Illinois project: recorded positions",
      subtitle = "Coordinate groups do not certify physical sites\nor independent samples",
      x = "Published WGS84 longitude", y = "Published WGS84 latitude") + ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(target, "recorded_positions.png"), map, width = 7, height = 5, dpi = 160)
  displayed <- result$declustering[result$declustering$element %in% c("Ag", "Pb", "Zn", "Cu", "Hg"), ]
  displayed$cell_km <- displayed$cell_size / 1000
  comparison <- ggplot2::ggplot(displayed, ggplot2::aes(x = cell_km, colour = element, group = element)) +
    ggplot2::geom_line(ggplot2::aes(y = geomean)) + ggplot2::geom_point(ggplot2::aes(y = geomean)) +
    ggplot2::geom_line(ggplot2::aes(y = unweighted_geomean), linetype = "dashed") +
    ggplot2::scale_y_log10() + ggplot2::scale_x_continuous(breaks = c(5, 10, 20)) +
    ggplot2::labs(title = "Declustering sensitivity: selected elements",
      subtitle = "Solid: weighted; dashed: unweighted on the same usable subset",
      caption = "Descriptive fits; no spatial confidence limits or source attribution",
      x = "Cell width (km); 100 grid placements", y = "Geometric mean (ug/kg)", colour = "Element") +
    ggplot2::theme_minimal()
  if (nrow(displayed)) ggplot2::ggsave(file.path(target, "declustering_sensitivity.png"), comparison, width = 8, height = 5, dpi = 160)
  provenance <- list(package_version = as.character(utils::packageVersion("geoCensor")),
    project_name = project, state = "IL", input_paths = normalizePath(files), input_md5 = before,
    scanned_metadata_records = main$n_input, scanned_chemistry_records = chemistry$n_input,
    projection = projection, cell_sizes_m = c(5000, 10000, 20000), n_offsets_per_axis = 10,
    observation_rule = "All recorded species for selected metadata IDs; no chemistry selection",
    statistical_status = "Descriptive; weighted limits absent; spatial limits experimental")
  saveRDS(provenance, file.path(target, "provenance.rds"))
  writeLines(capture.output(sessionInfo()), file.path(target, "session_info.txt"))
  writeLines(c("USGS Upper Illinois project worked example", capture.output(print(counts)),
    "All source fields are read as text. Source-record ordinals exclude the header and blank lines.",
    "Original coordinates are retained; projected coordinates are appended, not corrected.",
    "Recorded positions do not certify physical sites, replicates or independence.",
    "The cohort mixes unresolved sampling supports and is not a geological background population.",
    "All-left/unavailable species remain in estimates.csv; no nondetect or missing-data substitution.",
    "Ordinary intervals assume the working model/independence; weighted intervals are absent.",
    "No validated spatial intervals, ore/anomaly classification, or contamination attribution.",
    "No collection date, historical exposure or dilution factor is inferred.",
    "Use a new output directory for another run; partial outputs from a write failure are retained."),
    file.path(target, "README.txt"))
  outputs <- list.files(target, full.names = TRUE)
  readr::write_csv(data.frame(file = basename(outputs), bytes = file.info(outputs)$size,
    md5 = unname(tools::md5sum(outputs))), file.path(target, "output_checksums.csv"))
  message("Illinois worked example saved to ", normalizePath(target))
  invisible(list(results = result, support = support, provenance = provenance))
}
run_illinois(commandArgs(trailingOnly = TRUE))
