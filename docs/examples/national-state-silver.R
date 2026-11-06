#!/usr/bin/env Rscript
# National/state silver example. Download native NGDB exports separately.
# This repository script uses geoCensor 1.2.0; it is not an installed command.

silver_scan <- function(path, required, select, chunk_size = 10000L) {
  types <- readr::cols(.default = readr::col_character())
  header <- readr::read_csv(path, n_max = 0, col_types = types,
    na = character(), trim_ws = FALSE, name_repair = "check_unique",
    show_col_types = FALSE, progress = FALSE)
  if (!all(required %in% names(header))) stop("Missing columns in ", basename(path),
    ": ", paste(setdiff(required, names(header)), collapse = ", "), call. = FALSE)
  pieces <- list(); ordinals <- list(); scanned <- 0
  callback <- readr::SideEffectChunkCallback$new(function(x, pos) {
    if (!identical(names(x), names(header)) || nrow(readr::problems(x)) ||
        !all(vapply(x, is.character, logical(1))) || pos != scanned + 1) {
      stop("CSV parsing/schema problem in ", basename(path), " near record ", scanned + 1, call. = FALSE)
    }
    keep <- select(x)
    if (!is.logical(keep) || length(keep) != nrow(x) || anyNA(keep)) stop("Invalid selection mask")
    ii <- which(keep)
    if (length(ii)) {
      pieces[[length(pieces) + 1L]] <<- as.data.frame(x[ii, , drop = FALSE])
      ordinals[[length(ordinals) + 1L]] <<- scanned + ii
    }
    scanned <<- scanned + nrow(x)
  })
  withCallingHandlers(readr::read_csv_chunked(path, callback, chunk_size = chunk_size,
    col_types = types, na = character(), trim_ws = FALSE,
    show_col_types = FALSE, progress = FALSE), warning = function(w) {
      stop("CSV parsing warning: ", conditionMessage(w), call. = FALSE)
    })
  list(data = if (length(pieces)) as.data.frame(dplyr::bind_rows(pieces)) else as.data.frame(header),
    ordinals = as.numeric(unlist(ordinals, use.names = FALSE)), scanned = scanned)
}

silver_numeric_coordinate <- function(x) {
  missing <- is.na(x) | trimws(x) %in% c("", "NA")
  value <- suppressWarnings(as.numeric(x))
  if (any(!missing & (!is.finite(value) | is.na(value)))) stop("Malformed coordinate text", call. = FALSE)
  value[missing] <- NA_real_
  value
}

silver_project <- function(metadata) {
  lon <- silver_numeric_coordinate(metadata$long_wgs84)
  lat <- silver_numeric_coordinate(metadata$lat_wgs84)
  bounded <- is.finite(lon) & is.finite(lat) & abs(lon) <= 180 & abs(lat) <= 90
  valid <- bounded & !(lon == 0 & lat == 0)
  metadata$state_analysis <- toupper(trimws(metadata$state))
  metadata$national_component <- ifelse(metadata$state_analysis == "AK", "ALASKA",
    ifelse(metadata$state_analysis == "HI", "HAWAII",
      ifelse(metadata$state_analysis %in% c(setdiff(state.abb, c("AK", "HI")), "DC"), "CONUS", NA_character_)))
  for (name in c("x_state_m", "y_state_m", "x_national_m", "y_national_m")) metadata[[name]] <- NA_real_
  definitions <- list()
  project <- function(ii, crs) {
    if (!length(ii)) return(matrix(numeric(), ncol = 2))
    p <- sf::st_as_sf(data.frame(lon = lon[ii], lat = lat[ii]), coords = c("lon", "lat"), crs = 4326)
    xy <- sf::st_coordinates(sf::st_transform(p, crs))
    if (any(!is.finite(xy))) stop("Nonfinite projected coordinate", call. = FALSE)
    xy
  }
  # Centers use ALL explicitly-US metadata, before Ag selection, matching the
  # source study. (0,0) is retained for center derivation but never located Ag.
  for (state in c(state.abb, "DC")) {
    center_rows <- which(metadata$state_analysis == state & bounded)
    if (!length(center_rows)) next
    crs <- sprintf("+proj=laea +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m +no_defs",
      stats::median(lat[center_rows]), stats::median(lon[center_rows]))
    ii <- which(metadata$state_analysis == state & valid)
    xy <- project(ii, crs)
    metadata$x_state_m[ii] <- xy[, 1]; metadata$y_state_m[ii] <- xy[, 2]
    definitions[[length(definitions) + 1L]] <- data.frame(domain = state, role = "state",
      n_center_metadata = length(center_rows), n_projected_metadata = length(ii), crs = crs,
      wkt = sf::st_crs(crs)$wkt)
  }
  national <- c(CONUS = "EPSG:5070", ALASKA = "EPSG:3338",
    HAWAII = "+proj=laea +lat_0=20.5 +lon_0=-157 +datum=WGS84 +units=m +no_defs")
  for (component in names(national)) {
    ii <- which(metadata$national_component == component & valid)
    xy <- project(ii, national[[component]])
    metadata$x_national_m[ii] <- xy[, 1]; metadata$y_national_m[ii] <- xy[, 2]
    definitions[[length(definitions) + 1L]] <- data.frame(domain = component, role = "national",
      n_center_metadata = NA_integer_, n_projected_metadata = length(ii), crs = national[[component]],
      wkt = sf::st_crs(national[[component]])$wkt)
  }
  list(metadata = metadata, projections = as.data.frame(dplyr::bind_rows(definitions)))
}

silver_read_medium <- function(main_path, chemistry_path, medium) {
  if (!medium %in% c("soil", "sediment")) stop("Unknown medium")
  seen <- character()
  main <- silver_scan(main_path,
    c("lab_id", "country", "state", "primary_class", "lat_wgs84", "long_wgs84"), function(x) {
      ids <- trimws(x$lab_id)
      if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(ids %in% seen)) {
        stop("Metadata IDs must be nonblank and unique, including across chunks", call. = FALSE)
      }
      seen <<- c(seen, ids)
      toupper(trimws(x$country)) %in% c("US", "USA", "UNITED STATES", "U.S.")
    })
  if (any(c("metadata_source_record", "chemistry_source_record", "lab_id_source",
    "state_analysis", "national_component", "x_state_m", "y_state_m", "x_national_m", "y_national_m") %in% names(main$data))) {
    stop("Supply native metadata without example-generated columns", call. = FALSE)
  }
  main$data$metadata_source_record <- main$ordinals
  ids <- trimws(main$data$lab_id)
  chem <- silver_scan(chemistry_path, c("lab_id", "species", "qvalue", "unit", "amethod"),
    function(x) x$species == "Ag" & trimws(x$lab_id) %in% ids)
  if (any(c("chemistry_source_record", "metadata_source_record", "lab_id_source") %in% names(chem$data))) {
    stop("Supply native chemistry without example-generated columns", call. = FALSE)
  }
  if (anyDuplicated(trimws(chem$data$lab_id))) stop("Duplicate Ag best-value sample keys; review without averaging", call. = FALSE)
  chem$data$chemistry_source_record <- chem$ordinals
  projection <- silver_project(main$data)
  selected <- which(ids %in% trimws(chem$data$lab_id))
  original_metadata <- main$data[selected, , drop = FALSE]
  if (any(tolower(trimws(original_metadata$primary_class)) != medium)) {
    stop("Selected primary_class does not match the declared medium", call. = FALSE)
  }
  projected <- projection$metadata[selected, , drop = FALSE]
  original_chemistry <- chem$data
  chem$data$lab_id_source <- chem$data$lab_id
  chem$data$lab_id <- trimws(chem$data$lab_id)
  parsed <- geoCensor::build_element(chem$data, "Ag")
  if (nrow(parsed) != nrow(chem$data)) stop("Ag parsing changed row count")
  joined <- as.data.frame(geoCensor::join_soil_main(parsed, projected,
    extra_cols = setdiff(names(projected), "lab_id")))
  if (nrow(joined) != nrow(chem$data) || anyNA(joined$metadata_source_record) ||
      !identical(joined$chemistry_source_record, chem$data$chemistry_source_record)) stop("Source linkage changed")
  list(data = joined, original_metadata = original_metadata, original_chemistry = original_chemistry,
    projections = projection$projections,
    metadata_counts = table(factor(projection$metadata$state_analysis, levels = c(state.abb, "DC"))),
    inventory = data.frame(metric = c("metadata_scanned", "explicit_us_metadata", "chemistry_scanned",
      "explicit_us_Ag", "Ag_outside_50_states_dc"),
      n = c(main$scanned, nrow(main$data), chem$scanned, nrow(chem$data),
        sum(!joined$state_analysis %in% c(state.abb, "DC")))))
}

# Public parsed fields define intervals; this mirrors the API's documented
# eligibility rule. Zero here is a left interval endpoint, not an imputed assay.
silver_intervals <- function(d) {
  left <- as.logical(d$left_cen); right <- as.logical(d$right_cen)
  left[is.na(left)] <- FALSE; right[is.na(right)] <- FALSE
  if (any(left & right)) stop("Contradictory censor flags")
  il <- left & is.finite(d$left_limit) & d$left_limit > 0
  ir <- right & is.finite(d$right_limit) & d$right_limit > 0
  ie <- !left & !right & is.finite(d$detected_val) & d$detected_val > 0
  lo <- hi <- rep(NA_real_, nrow(d))
  lo[il] <- 0; hi[il] <- d$left_limit[il]
  lo[ir] <- d$right_limit[ir]; hi[ir] <- Inf
  lo[ie] <- hi[ie] <- d$detected_val[ie]
  list(ok = il | ir | ie, lower = lo, upper = hi)
}

# At every placement each occupied cell across ALL projection components gets
# equal mass. Component IDs prevent cells in Alaska and CONUS from colliding.
# National anchor=fixed0 reproduces the scale study. State anchor=minimum
# reproduces the package cell_declustering_weights() convention.
silver_weights <- function(d, size_m, n_offsets = 10L, anchor = c("fixed0", "minimum")) {
  anchor <- match.arg(anchor)
  if (!is.numeric(size_m) || length(size_m) != 1L || !is.finite(size_m) || size_m <= 0 ||
      length(n_offsets) != 1L || !is.finite(n_offsets) || n_offsets < 1 || n_offsets != as.integer(n_offsets)) stop("Invalid grid configuration")
  if (!all(c("component", "x_m", "y_m") %in% names(d)) ||
      anyNA(d$component) || any(!nzchar(d$component)) || any(!is.finite(d$x_m)) || any(!is.finite(d$y_m))) stop("Weights require located records and components")
  if (!nrow(d)) return(list(weights = numeric(), cells_min = 0L, cells_mean = 0, cells_max = 0L))
  groups <- split(seq_len(nrow(d)), d$component, drop = TRUE)
  xy <- lapply(groups, function(ii) {
    x <- d$x_m[ii]; y <- d$y_m[ii]
    if (anchor == "minimum") { x <- x - min(x); y <- y - min(y) }
    cbind(x / size_m, y / size_m)
  })
  w <- numeric(nrow(d)); occupied <- integer(); offset <- (seq_len(n_offsets) - 1) / n_offsets
  for (dx in offset) for (dy in offset) {
    indexes <- counts <- list(); total <- 0L
    for (component in names(groups)) {
      ii <- floor(xy[[component]][, 1] + dx)
      jj <- floor(xy[[component]][, 2] + dy)
      # Signed cell pairs are encoded as strings without integer collisions.
      key <- paste(sprintf("%.0f", ii), sprintf("%.0f", jj), sep = ":")
      index <- match(key, unique(key)); count <- tabulate(index)
      indexes[[component]] <- index; counts[[component]] <- count
      total <- total + length(count)
    }
    occupied <- c(occupied, total)
    for (component in names(groups)) w[groups[[component]]] <- w[groups[[component]]] +
      1 / (total * counts[[component]][indexes[[component]]])
  }
  w <- w / sum(w)
  stopifnot(all(is.finite(w)), all(w > 0), abs(sum(w) - 1) < 1e-12)
  list(weights = w, cells_min = min(occupied), cells_mean = mean(occupied), cells_max = max(occupied))
}

silver_fit <- function(d, weights = NULL) {
  if (!is.null(weights)) {
    if (length(weights) != nrow(d) || any(!is.finite(weights)) || any(weights < 0)) stop("Invalid weights")
    d$.silver_weight <- weights
  }
  fit <- as.data.frame(geoCensor::flag_estimate_information(geoCensor::fit_cen_lognorm(d,
    weight_col = if (is.null(weights)) NULL else ".silver_weight")))
  iv <- silver_intervals(d)
  w <- if (is.null(weights)) rep(1, nrow(d)) else weights
  use <- iv$ok & w > 0
  if (fit$n_total != sum(use)) stop("Example/API interval eligibility disagreement")
  w <- w[use]; if (length(w)) w <- w / sum(w)
  lo <- iv$lower[use]; hi <- iv$upper[use]
  # Point estimates only: even the ordinary intervals do not describe spatial
  # or preferential-sampling uncertainty, so omit every confidence-limit field.
  fit$ci_low <- fit$ci_high <- NA_real_
  fit$geomean_mgkg <- fit$geomean / 1000
  fit$meanlog_mgkg <- fit$sdlog <- fit$arithmetic_raw_mgkg <- fit$arithmetic_mean_mgkg <- NA_real_
  fit$arithmetic_status <- "estimate_unavailable"
  fit$scale_status <- "estimate_unavailable"
  if (is.finite(fit$geomean) && fit$geomean > 0) {
    # The public API returns GM; recover scale using the same likelihood and
    # confirm its location against the public result before exporting moments.
    native <- tryCatch(withCallingHandlers(survival::survreg(
      survival::Surv(lo, hi, type = "interval2") ~ 1, dist = "lognormal", weights = w * length(w)),
      warning = function(e) stop(conditionMessage(e))), error = function(e) e)
    if (inherits(native, "error")) {
      fit$scale_status <- paste0("scale_refit_failed: ", conditionMessage(native))
    } else {
      mu <- unname(stats::coef(native)[1]); sigma <- native$scale
      if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0 || abs(log(fit$geomean) - mu) > 1e-7) {
        stop("Scale refit disagrees with the public geometric mean", call. = FALSE)
      }
      fit$scale_status <- "same_likelihood_checked"
      fit$meanlog_mgkg <- mu - log(1000); fit$sdlog <- sigma
      fit$arithmetic_raw_mgkg <- exp(mu + sigma^2 / 2) / 1000
      fit$arithmetic_status <- if (!is.finite(fit$arithmetic_raw_mgkg)) "nonfinite_model_moment" else
        if (fit$arithmetic_raw_mgkg > 1e6) "model_mean_exceeds_100_percent_mass" else "model_implied_lognormal_mean"
      if (fit$arithmetic_status == "model_implied_lognormal_mean") fit$arithmetic_mean_mgkg <- fit$arithmetic_raw_mgkg
    }
  }
  fit$sample_arithmetic_lower_mgkg <- if (length(w)) sum(w * lo) / 1000 else NA_real_
  fit$sample_arithmetic_upper_mgkg <- if (length(w)) sum(w * hi) / 1000 else NA_real_
  fit$model_mean_below_sample_lower <- if (is.finite(fit$arithmetic_raw_mgkg))
    fit$arithmetic_raw_mgkg < fit$sample_arithmetic_lower_mgkg else NA
  fit$weighted_exact_fraction <- if (length(w)) sum(w[lo == hi]) else NA_real_
  fit$weighted_left_fraction <- if (length(w)) sum(w[lo == 0]) else NA_real_
  fit$weighted_right_fraction <- if (length(w)) sum(w[is.infinite(hi)]) else NA_real_
  fit$weight_effective_n <- if (length(w)) 1 / sum(w^2) else 0
  fit$n_method_codes <- length(unique(d$amethod[use]))
  fit$method_codes <- paste(sort(unique(d$amethod[use])), collapse = ";")
  fit$spatial_confidence_status <- "not_calculated_experimental"
  fit
}

silver_analyze <- function(input, medium, directory, national_sizes = c(20000, 40000, 80000),
    state_sizes = c(5000, 10000, 20000), n_offsets = 10L) {
  d <- input$data; fits <- list(); case_id <- 0L
  weight_path <- file.path(directory, "row_weights.csv")
  utils::write.csv(data.frame(case_id = character(), chemistry_source_record = numeric(),
    cell_size_m = numeric(), weight = numeric()), weight_path, row.names = FALSE)
  add <- function(rows, domain, domain_type) {
    methods <- sort(unique(d$amethod[rows]))
    variants <- c(list(rows), lapply(methods, function(method) rows[d$amethod[rows] == method]))
    for (v in seq_along(variants)) {
      bucket <- d[variants[[v]], , drop = FALSE]; case_id <<- case_id + 1L
      id <- sprintf("%s_%04d", medium, case_id)
      bucket$x_m <- if (domain_type == "state") bucket$x_state_m else bucket$x_national_m
      bucket$y_m <- if (domain_type == "state") bucket$y_state_m else bucket$y_national_m
      bucket$component <- if (domain_type == "state") rep(domain, nrow(bucket)) else bucket$national_component
      valid <- silver_intervals(bucket)$ok
      located <- valid & is.finite(bucket$x_m) & is.finite(bucket$y_m) & !is.na(bucket$component)
      spatial <- bucket[located, , drop = FALSE]
      wanted <- if (domain_type == "state") domain else if (domain == "US50_DC") c(state.abb, "DC") else c(setdiff(state.abb, c("AK", "HI")), "DC")
      prefix <- data.frame(case_id = id, medium = medium, domain_type = domain_type, domain = domain,
        analysis_role = if (v == 1L) "ALL_METHODS_REFERENCE" else "METHOD_STRATUM",
        method_code = if (v == 1L) NA_character_ else methods[v - 1L],
        n_metadata = sum(input$metadata_counts[wanted]),
        n_input = nrow(bucket), n_usable_chemistry = sum(valid), n_located = sum(located),
        n_unlocated_usable = sum(valid & !located),
        n_coordinate_pairs = nrow(unique(spatial[, c("component", "x_m", "y_m"), drop = FALSE])),
        data_status = if (!nrow(bucket)) "no_reported_Ag" else if (!any(located)) "no_located_usable_Ag" else "observed_Ag_records")
      save_fit <- function(fit, basis, size = NA_real_, grid = NULL) {
        fits[[length(fits) + 1L]] <<- cbind(prefix, estimate_basis = basis, cell_size_m = size,
          n_origins = if (is.null(grid)) 0L else n_offsets^2,
          cells_min = if (is.null(grid)) NA_real_ else grid$cells_min,
          cells_mean = if (is.null(grid)) NA_real_ else grid$cells_mean,
          cells_max = if (is.null(grid)) NA_real_ else grid$cells_max, fit)
      }
      save_fit(silver_fit(bucket), "unweighted_all")
      save_fit(silver_fit(spatial), "unweighted_located")
      sizes <- if (domain_type == "state") state_sizes else national_sizes
      for (size in sizes) {
        grid <- silver_weights(spatial, size, n_offsets,
          anchor = if (domain_type == "state") "minimum" else "fixed0")
        save_fit(silver_fit(spatial, grid$weights), "cell_declustered", size, grid)
        if (nrow(spatial)) utils::write.table(data.frame(case_id = id,
          chemistry_source_record = spatial$chemistry_source_record, cell_size_m = size, weight = grid$weights),
          weight_path, append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",", qmethod = "double")
      }
      message(medium, ": ", domain, " / ", if (v == 1L) "pooled reference" else methods[v - 1L],
        " (", nrow(spatial), " located usable records)")
    }
  }
  for (state in c(state.abb, "DC")) add(which(d$state_analysis == state), state, "state")
  add(which(d$state_analysis %in% c(state.abb, "DC")), "US50_DC", "national")
  add(which(d$state_analysis %in% c(setdiff(state.abb, c("AK", "HI")), "DC")), "CONUS", "national")
  estimates <- as.data.frame(dplyr::bind_rows(fits))
  write <- function(x, name) utils::write.csv(x, file.path(directory, name), row.names = FALSE, na = "NA")
  write(estimates, "estimates.csv")
  write(estimates[estimates$domain_type == "state" & estimates$analysis_role == "ALL_METHODS_REFERENCE", ], "state_summary.csv")
  write(estimates[estimates$domain_type == "national" & estimates$analysis_role == "ALL_METHODS_REFERENCE", ], "national_summary.csv")
  saveRDS(list(estimates = estimates, configuration = list(national_sizes_m = national_sizes,
    state_sizes_m = state_sizes, offsets_per_axis = n_offsets)), file.path(directory, "results.rds"))
  invisible(estimates)
}

silver_run <- function(args = commandArgs(TRUE)) {
  usage <- paste("Usage: Rscript docs/examples/national-state-silver.R",
    "soil_main.csv soil_bestvalue.csv sediment_main.csv sediment_bestvalue.csv NEW_output_directory")
  if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  if (length(args) != 5L || any(startsWith(args, "--"))) stop(usage, call. = FALSE)
  for (pkg in c("geoCensor", "sf")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " first", call. = FALSE)
  if (as.character(utils::packageVersion("geoCensor")) != "1.2.0") stop("This example requires geoCensor 1.2.0", call. = FALSE)
  files <- stats::setNames(args[1:4], c("soil_main", "soil_bestvalue", "sediment_main", "sediment_bestvalue"))
  if (any(!file.exists(files) | dir.exists(files))) stop("Inputs must be readable native CSV files", call. = FALSE)
  files <- stats::setNames(normalizePath(files, mustWork = TRUE), names(files))
  if (anyDuplicated(files)) stop("Supply four distinct input files", call. = FALSE)
  target <- args[5]
  check_target <- function() {
    link <- Sys.readlink(target)
    if (file.exists(target) || (!is.na(link) && nzchar(link))) stop("Output already exists; choose a new directory", call. = FALSE)
    if (!dir.exists(dirname(target))) stop("Output parent directory must exist", call. = FALSE)
  }
  check_target()
  before <- tools::md5sum(files)
  inputs <- list(
    soil = silver_read_medium(files[1], files[2], "soil"),
    sediment = silver_read_medium(files[3], files[4], "sediment"))
  if (!identical(unname(before), unname(tools::md5sum(files)))) stop("Inputs changed during reading", call. = FALSE)
  check_target()
  if (!dir.create(target)) stop("Cannot create output directory", call. = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(target, "sessionInfo.txt"))
  writeLines(capture.output(sf::sf_extSoftVersion()), file.path(target, "projection_software.txt"))
  utils::write.csv(data.frame(input = names(files), file = unname(files), bytes = file.info(files)$size,
    md5 = unname(before)), file.path(target, "input_checksums.csv"), row.names = FALSE)
  configuration <- list(package = "geoCensor", package_version = "1.2.0", element = "Ag",
    national_domains = c("US50_DC", "CONUS"), national_cell_sizes_m = c(20000, 40000, 80000),
    state_cell_sizes_m = c(5000, 10000, 20000), offsets_per_axis = 10L, n_origins = 100L,
    national_anchor = "fixed0 in each projection component", state_anchor = "per-case component minimum extent",
    country_labels = c("US", "USA", "UNITED STATES", "U.S."),
    state_projection_center = "coordinate medians of all bounded explicitly-US metadata in each state, rounded to 10 decimals",
    target = "occupied sampled cells; pooled method reference and separate observed method strata",
    uncertainty = "No spatial confidence intervals; cell widths and origins are sensitivity settings")
  saveRDS(configuration, file.path(target, "configuration.rds"))
  writeLines(capture.output(dput(configuration)), file.path(target, "configuration.txt"))
  for (medium in names(inputs)) {
    out <- file.path(target, medium); dir.create(out)
    input <- inputs[[medium]]
    for (name in c("original_metadata", "original_chemistry", "projections", "inventory")) {
      utils::write.csv(input[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "NA")
    }
    saveRDS(input$data, file.path(out, "analysis.rds"))
    support_fields <- intersect(c("primary_class", "sample_source", "method_collected", "mesh_size", "horizon", "depth", "date_collect", "prep"), names(input$data))
    support <- dplyr::bind_rows(lapply(support_fields, function(field) {
      z <- table(input$data[[field]], useNA = "ifany")
      data.frame(field = field, value = names(z), n = as.numeric(z))
    }))
    utils::write.csv(support, file.path(out, "sampling_support_inventory.csv"), row.names = FALSE, na = "NA")
    silver_analyze(input, medium, out)
  }
  if (!identical(unname(before), unname(tools::md5sum(files)))) stop("Inputs changed during analysis; inspect retained partial output", call. = FALSE)
  artifacts <- list.files(target, full.names = TRUE, recursive = TRUE)
  utils::write.csv(data.frame(file = substring(artifacts, nchar(target) + 2L), bytes = file.info(artifacts)$size,
    md5 = unname(tools::md5sum(artifacts))), file.path(target, "output_checksums.csv"), row.names = FALSE)
  writeLines("Completed example computation. Computational completion is not scientific validation.", file.path(target, "COMPLETED.txt"))
  invisible(target)
}

if (sys.nframe() == 0L) silver_run()
