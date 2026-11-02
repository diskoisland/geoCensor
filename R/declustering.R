#' Calculate two-dimensional cell-declustering weights
#'
#' Each occupied square cell receives equal total weight, shared among its
#' observations. Weights are averaged over a deterministic lattice of grid
#' origins to reduce dependence on grid placement.
#'
#' @param df Data frame of observations from one spatial domain.
#' @param x_col,y_col Names of numeric projected coordinate columns. Longitude
#'   and latitude in degrees must be projected before use.
#' @param cell_size Positive square-cell width in the coordinate units.
#' @param n_offsets Number of equally spaced origins per axis. The default
#'   uses 100 origins (10 by 10); no random numbers are generated.
#' @return Numeric weights aligned with input rows, summing to one when valid
#'   coordinates exist. Nonfinite coordinates receive zero weight. If none
#'   exist, all weights are zero.
#' @details
#' Apply separately within a scientifically defined domain and sampling medium.
#' Coordinates must use an appropriate projected coordinate system; this
#' function cannot infer or verify the coordinate reference system from numbers.
#' Repeated locations share influence; observations are not averaged or deleted.
#' Weights describe the occupied cells and do not represent unsampled terrain.
#' Cell size should be assessed against sampling spacing and geology rather
#' than automatically selected to minimize the resulting concentration estimate.
#' @references
#' Deutsch, C. V. (2015). Cell Declustering Parameter Selection.
#' \url{https://geostatisticslessons.com/lessons/celldeclustering}
#' @export
cell_declustering_weights <- function(df, x_col, y_col, cell_size, n_offsets = 10) {
  if (!all(c(x_col, y_col) %in% names(df))) stop("Missing projected coordinate columns")
  x <- df[[x_col]]
  y <- df[[y_col]]
  if (!is.numeric(x) || !is.numeric(y)) stop("Projected coordinates must be numeric")
  if (!is.numeric(cell_size) || length(cell_size) != 1L ||
      !is.finite(cell_size) || cell_size <= 0) stop("cell_size must be one positive number")
  if (!is.numeric(n_offsets) || length(n_offsets) != 1L ||
      !is.finite(n_offsets) || n_offsets < 1 || n_offsets != floor(n_offsets)) {
    stop("n_offsets must be one positive whole number")
  }
  valid <- is.finite(x) & is.finite(y)
  weights <- numeric(nrow(df))
  if (!any(valid)) return(weights)
  # Anchor to the data extent so changing coordinate origin does not change weights.
  x <- (x[valid] - min(x[valid])) / cell_size
  y <- (y[valid] - min(y[valid])) / cell_size
  if (any(!is.finite(x)) || any(!is.finite(y)) || any(x > 2^52) || any(y > 2^52)) {
    stop("cell_size is too small relative to the coordinate extent")
  }
  averaged <- numeric(length(x))
  offsets <- (seq_len(n_offsets) - 1) / n_offsets
  for (dx in offsets) {
    for (dy in offsets) {
      ix <- floor(x + dx)
      iy <- floor(y + dy)
      width <- max(ix) + 1
      # Numeric cell keys avoid allocating strings for large national datasets.
      # Fall back to strings when the combined key would lose integer precision.
      cells <- if (width * (max(iy) + 1) <= 2^52) ix + width * iy else paste(ix, iy, sep = ":")
      index <- match(cells, unique(cells))
      counts <- tabulate(index)
      averaged <- averaged + 1 / (length(counts) * counts[index])
    }
  }
  weights[valid] <- averaged / sum(averaged)
  weights
}

#' Assess censored lognormal estimates across declustering cell sizes
#'
#' Compares a common unweighted baseline with spatially weighted censored
#' lognormal estimates. Weights use only rows with usable chemistry and finite
#' projected coordinates, so the comparison uses the same observations.
#'
#' @param df Data frame containing the standard censoring columns used by
#'   [fit_cen_lognorm()] and projected coordinates.
#' @param cell_sizes Positive cell widths to assess, in coordinate units.
#' @param x_col,y_col Names of numeric projected coordinate columns.
#' @param n_offsets Grid origins per axis, passed to [cell_declustering_weights()].
#' @return A tibble with one row per cell size, original and usable row counts,
#'   the common unweighted geometric mean, weighted fit fields, exact_fraction,
#'   weighted_exact_fraction, weighted_left_fraction, weighted_right_fraction,
#'   and Kish weight effective sample size. The latter measures weight imbalance, not the number
#'   of spatially independent observations. Confidence limits are missing because
#'   a spatial uncertainty procedure has not been applied.
#' @details
#' This is a sensitivity analysis, not a complete correction for preferential
#' sampling. It does not infer concentrations in unsampled cells. Apply separately
#' within each sampling medium and meaningful spatial/geological domain.
#' exact_fraction is the proportion of usable rows with exact observations.
#' The weighted fractions describe the share of total declustering weight carried
#' by each observation type on that same subset; they sum to one when data exist
#' and are missing when no rows are usable. They are descriptive, not effective
#' information counts or validated reliability thresholds. A small share of weight
#' on exact values can make location and scale estimation sensitive even when
#' many raw records are present. No new exclusion rule or confidence interval
#' follows from these fractions.
#' @export
declustering_sensitivity <- function(df, cell_sizes, x_col, y_col, n_offsets = 10) {
  if (!is.numeric(cell_sizes) || !length(cell_sizes) ||
      any(!is.finite(cell_sizes)) || any(cell_sizes <= 0)) {
    stop("cell_sizes must contain positive finite numbers")
  }
  intervals <- .censored_intervals(
    df, "detected_val", "left_cen", "left_limit", "right_cen", "right_limit",
    "declustering_sensitivity"
  )
  # Validate the coordinate arguments even when there are no usable rows.
  spatial_weights <- cell_declustering_weights(df, x_col, y_col, cell_sizes[1], n_offsets)
  usable <- intervals$ok & spatial_weights > 0
  subset <- df[usable, , drop = FALSE]
  exact <- intervals$lower[usable] == intervals$upper[usable]
  left <- intervals$lower[usable] == 0
  right <- is.infinite(intervals$upper[usable])
  baseline <- fit_cen_lognorm(subset)
  dplyr::bind_rows(lapply(cell_sizes, function(size) {
    weights <- cell_declustering_weights(subset, x_col, y_col, size, n_offsets)
    weighted_data <- subset
    weighted_data$.declustering_weight <- weights
    fit <- fit_cen_lognorm(weighted_data, weight_col = ".declustering_weight")
    fraction <- function(index) {
      if (length(weights) && sum(weights) > 0) sum(weights[index]) / sum(weights) else NA_real_
    }
    tibble::tibble(
      cell_size = size, n_origins = n_offsets^2,
      n_input = nrow(df), n_excluded = sum(!usable),
      unweighted_geomean = baseline$geomean,
      exact_fraction = if (length(exact)) mean(exact) else NA_real_,
      weighted_exact_fraction = fraction(exact),
      weighted_left_fraction = fraction(left),
      weighted_right_fraction = fraction(right),
      weight_effective_n = if (length(weights) && sum(weights) > 0) 1 / sum(weights^2) else 0
    ) %>% dplyr::bind_cols(fit)
  }))
}
