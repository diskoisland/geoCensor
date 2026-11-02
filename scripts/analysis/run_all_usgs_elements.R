#!/usr/bin/env Rscript
# Reproducible full-data run. Uses local package source and installed sf/data.table.
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  pkgload::load_all('.', quiet = TRUE)
})
setDTthreads(2)
options(stringsAsFactors = FALSE)
output_root <- file.path('output', 'usgs_all_elements_2026-10-31')
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
cell_sizes <- c(5000, 10000, 20000)
writeLines(capture.output(sessionInfo()), file.path(output_root, 'session_info.txt'))
writeLines(c(
  'Scope: samples explicitly identified as United States in each NGDB main.csv.',
  'Soil and sediment are analyzed separately; all reported species are retained.',
  'Pooled estimates are unweighted. Grouped analyses use reported state/territory.',
  'Unknown states are retained as UNKNOWN_STATE; no coordinates are assigned to that group.',
  'Each state/territory uses its own local WGS84 Lambert azimuthal equal-area projection.',
  'Cell sizes: 5000, 10000, 20000 metres; 100 deterministic grid origins per size.',
  'These exploratory scales are not selected as optimal geological sampling spacings.',
  'Spatial comparisons use matched usable chemistry/coordinate subsets within each element/group.',
  'Distinct lab_id samples at the same location remain separate and share cell influence.',
  'Unsupported units (e.g. grams rather than concentration) remain excluded and counted.',
  'Rounded upper-censoring suffixes are recovered from the first summary value when it agrees with rounded qvalue.',
  'Different horizons and size fractions are retained as reported; these are not design-unbiased population estimates.',
  'Weighted confidence limits are omitted; weight effective n does not measure spatial independence.',
  'Source CSVs are read only. Results are checkpointed by medium and species.'
), file.path(output_root, 'analysis_notes.txt'))
for (medium in c('soil', 'sed')) {
  cat('\nREADING', medium, '\n'); flush.console()
  root <- file.path('..', 'data', paste0('ngdb', medium))
  input <- file.path(root, 'bestvalue.csv')
  chemistry <- fread(input, select = c('lab_id','species','unit','qvalue','summary'),
                      colClasses = c(lab_id = 'character', qvalue = 'character'))
  metadata <- fread(file.path(root, 'main.csv'),
                    select = c('lab_id','state','country','lat_wgs84','long_wgs84'),
                    colClasses = c(lab_id = 'character'))
  input_n <- nrow(chemistry)
  meta_n <- nrow(metadata)
  if (anyDuplicated(metadata$lab_id)) stop('Duplicate metadata IDs')
  metadata <- metadata[toupper(trimws(country)) %in% c('US','USA','UNITED STATES','U.S.')]
  metadata[, state := toupper(trimws(state))]
  metadata[is.na(state) | state == '', state := 'UNKNOWN_STATE']
  chemistry <- chemistry[lab_id %in% metadata$lab_id]
  if (anyDuplicated(chemistry[, .(lab_id, species)])) stop('Duplicate sample/species keys')
  folder <- file.path(output_root, medium)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(folder, 'checkpoints'), showWarnings = FALSE)
  fwrite(data.table(medium = medium, input_chemistry_rows = input_n,
                     us_chemistry_rows = nrow(chemistry), excluded_non_us_or_unmatched_rows = input_n - nrow(chemistry),
                     input_metadata_rows = meta_n, us_metadata_rows = nrow(metadata),
                     n_species = uniqueN(chemistry$species), input_path = normalizePath(input)),
          file.path(folder, 'input_inventory.csv'))
  metadata[, c('x_m','y_m') := .(NA_real_, NA_real_)]
  projections <- list()
  for (region in unique(metadata$state)) {
    rows <- which(metadata$state == region & is.finite(metadata$lat_wgs84) &
                    is.finite(metadata$long_wgs84) & abs(metadata$lat_wgs84) <= 90 &
                    abs(metadata$long_wgs84) <= 180)
    if (!length(rows) || region == 'UNKNOWN_STATE') next
    lat0 <- median(metadata$lat_wgs84[rows])
    lon0 <- median(metadata$long_wgs84[rows])
    crs <- sprintf('+proj=laea +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m +no_defs', lat0, lon0)
    points <- st_as_sf(as.data.frame(metadata[rows]), coords = c('long_wgs84','lat_wgs84'), crs = 4326)
    xy <- st_coordinates(st_transform(points, crs))
    set(metadata, rows, 'x_m', xy[,1]); set(metadata, rows, 'y_m', xy[,2])
    projections[[region]] <- data.table(state = region, n_coordinates = length(rows), crs = crs)
  }
  fwrite(rbindlist(projections), file.path(folder, 'projection_definitions.csv'))
  fwrite(metadata, file.path(folder, 'metadata_projected.csv'))
  elements <- sort(unique(chemistry$species))
  for (i in seq_along(elements)) {
    element <- elements[i]
    checkpoint <- file.path(folder, 'checkpoints', sprintf('%03d.rds', i))
    if (file.exists(checkpoint)) {
      cached <- readRDS(checkpoint)
      if (!identical(cached$element, element)) stop('Checkpoint species mismatch')
      cat(medium, i, '/', length(elements), element, 'cached\n'); flush.console(); next
    }
    cat(medium, i, '/', length(elements), element, 'starting\n'); flush.console()
    subset <- as.data.frame(chemistry[species == element])
    start <- proc.time()[3]
    warnings <- character()
    result <- withCallingHandlers({
      list(
        pooled = analyze_elements(subset),
        state = analyze_elements(subset, main_df = as.data.frame(metadata), group_cols = 'state',
                                  cell_sizes = cell_sizes, x_col = 'x_m', y_col = 'y_m')
      )
    }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')})
    saveRDS(list(element = element, result = result, warnings = unique(warnings),
                 elapsed_seconds = unname(proc.time()[3] - start),
                 n_upper_limits_recovered = sum(geoCensor:::.qvalue_for_censoring(subset)$source == 'summary_upper_limit')), checkpoint)
    cat(medium, i, '/', length(elements), element, 'complete in', round(proc.time()[3] - start, 1), 'seconds\n'); flush.console()
  }
  checkpoints <- lapply(seq_along(elements), function(i) readRDS(file.path(folder, 'checkpoints', sprintf('%03d.rds', i))))
  for (scope in c('pooled','state')) {
    for (table_name in c('estimates','distributions','validation','declustering','joins')) {
      table <- dplyr::bind_rows(lapply(checkpoints, function(x) x$result[[scope]][[table_name]]))
      if (ncol(table)) fwrite(as.data.table(table), file.path(folder, paste0(scope,'_',table_name,'.csv')))
    }
  }
  fwrite(data.table(element = elements, elapsed_seconds = vapply(checkpoints, function(x) x$elapsed_seconds, numeric(1)),
                     n_upper_limits_recovered = vapply(checkpoints, function(x) as.numeric(x$n_upper_limits_recovered), numeric(1)),
                     warnings = vapply(checkpoints, function(x) paste(x$warnings, collapse = ' | '), character(1))),
          file.path(folder,'run_status.csv'))
  saveRDS(lapply(checkpoints, function(x) x$result), file.path(folder,'analysis.rds'))
  cat('FINISHED', medium, '\n'); flush.console()
  rm(chemistry, metadata, checkpoints); gc()
}
writeLines('Completed both soil and sediment analyses.', file.path(output_root,'COMPLETED.txt'))
