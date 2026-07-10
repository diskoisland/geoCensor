#' Summarise spatial coverage by state
#'
#' Counts valid coordinate records and occupied latitude/longitude grid cells by
#' state. This is a lightweight diagnostic for whether a state-level estimate is
#' supported by broad spatial coverage or by many samples clustered in fewer
#' places.
#'
#' @param df Data frame with state and coordinate columns.
#' @param state_col Name of the state column. Defaults to `"state"`.
#' @param lat_col Name of the latitude column. Defaults to `"lat_wgs84"`.
#' @param lon_col Name of the longitude column. Defaults to `"long_wgs84"`.
#' @param grid_deg Grid-cell size in decimal degrees. Defaults to `0.5`.
#'
#' @return A tibble with one row per state and columns for valid coordinate
#'   count, occupied grid-cell count, average samples per occupied cell, latitude
#'   and longitude ranges, and bounding-box area in square degrees.
#' @details
#' `spatial_coverage_by_state()` is a screening metric, not an area-weighted
#' sampling-design estimator. Larger `occupied_grid_cells` and lower
#' `samples_per_grid_cell` generally indicate broader spatial spread for a given
#' sample count. The `bbox_area_deg2` column is useful for spotting very narrow
#' spatial footprints, but it should not be interpreted as land area.
#' @export
spatial_coverage_by_state <- function(df,
                                      state_col = "state",
                                      lat_col = "lat_wgs84",
                                      lon_col = "long_wgs84",
                                      grid_deg = 0.5) {
  req <- c(state_col, lat_col, lon_col)
  missing <- req[!req %in% names(df)]
  if (length(missing) > 0) {
    stop(
      "spatial_coverage_by_state: Missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  if (!is.numeric(grid_deg) || length(grid_deg) != 1 ||
      !is.finite(grid_deg) || grid_deg <= 0) {
    stop("spatial_coverage_by_state: grid_deg must be one positive number")
  }

  empty <- tibble::tibble(
    state = character(),
    n_spatial = integer(),
    grid_deg = double(),
    occupied_grid_cells = integer(),
    samples_per_grid_cell = double(),
    lat_min = double(),
    lat_max = double(),
    lon_min = double(),
    lon_max = double(),
    lat_range = double(),
    lon_range = double(),
    bbox_area_deg2 = double()
  )

  points <- tibble::tibble(
    state = as.character(df[[state_col]]),
    lat = as.numeric(df[[lat_col]]),
    lon = as.numeric(df[[lon_col]])
  ) %>%
    dplyr::filter(
      !is.na(.data$state),
      is.finite(.data$lat),
      is.finite(.data$lon),
      .data$lat >= -90,
      .data$lat <= 90,
      .data$lon >= -180,
      .data$lon <= 180
    )

  if (nrow(points) == 0) {
    return(empty)
  }

  points %>%
    dplyr::mutate(
      lat_cell = floor(.data$lat / grid_deg) * grid_deg,
      lon_cell = floor(.data$lon / grid_deg) * grid_deg,
      grid_cell = paste(.data$lat_cell, .data$lon_cell, sep = "_")
    ) %>%
    dplyr::group_by(.data$state) %>%
    dplyr::summarise(
      n_spatial = dplyr::n(),
      occupied_grid_cells = dplyr::n_distinct(.data$grid_cell),
      samples_per_grid_cell = .data$n_spatial / .data$occupied_grid_cells,
      lat_min = min(.data$lat),
      lat_max = max(.data$lat),
      lon_min = min(.data$lon),
      lon_max = max(.data$lon),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      grid_deg = grid_deg,
      .after = "n_spatial"
    ) %>%
    dplyr::mutate(
      lat_range = .data$lat_max - .data$lat_min,
      lon_range = .data$lon_max - .data$lon_min,
      bbox_area_deg2 = .data$lat_range * .data$lon_range
    ) %>%
    dplyr::arrange(.data$state)
}
