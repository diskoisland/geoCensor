#' Audit geographic overlap of recorded analytical methods
#'
#' Counts occupied projected grid cells by method and shared cells between
#' named method codes at several cell sizes and deterministic grid placements.
#' Use one element and a scientifically selected sampling medium/domain.
#'
#' @param df Observations with the standard censoring columns used by
#'   [fit_cen_lognorm()] and projected coordinates. If species is supplied, its
#'   nonblank values must identify at most one element; select an element first.
#' @param cell_sizes Positive finite square-cell widths in coordinate units.
#'   Sorted unique sizes are returned.
#' @param x_col,y_col Distinct numeric projected coordinate column names.
#'   Project longitude/latitude first; a CRS cannot be verified from numbers.
#' @param n_offsets Equally spaced grid placements per axis. Default 2 evaluates
#'   four placements. Increase to inspect more placements; this is not a random
#'   sample or a confidence calculation.
#' @param method_col Method-code column; default amethod. An absent column leaves
#'   method information missing, rather than inferred.
#' @return A list with support (one row of input/chemistry/coordinate counts and
#'   the common grid anchor), methods (counts and occupied-cell minima/means/
#'   maxima per code/size), and pairs (shared/union-cell and Jaccard summaries
#'   across placements for each pair of named codes/size). Jaccard is shared
#'   occupied cells divided by their union. Counts describe records and exact
#'   coordinate pairs, not independently sampled sites.
#' @details
#' Uses only interval-valid chemistry with finite coordinate pairs. Missing
#' coordinates and excluded chemistry are counted separately. Every recorded
#' method bucket is retained in methods, including invalid-only or unlocated
#' buckets. Case and surrounding whitespace are normalized as in
#' [flag_analytical_methods()]. Blank/absent methods and literal UNKNOWN remain
#' explicit unresolved buckets and do not enter named-method pairs. Other
#' supplied codes are labels, not verified quantitative procedures.
#'
#' Pairs are retained even if one named code has no usable spatial observations;
#' their comparison metrics stay NA with zero evaluated placements. A code with
#' no usable spatial observations has missing occupied-cell metrics, not zero.
#' When both codes have spatial data, an observed zero shared-cell count is a
#' genuine zero for that placement/scale. Empty domains keep missing anchors
#' and empty method/pair tables. No location, concentration or method is imputed.
#'
#' All methods use the same grid at each placement, anchored to the minima of
#' all spatially usable observations, including unresolved methods. Cells are
#' half-open: floor((coordinate - anchor) / cell_size + offset), with offsets
#' 0, 1/n_offsets, ..., (n_offsets-1)/n_offsets. Anchoring makes the calculation
#' invariant to a common coordinate translation within numerical precision.
#' Minima/maxima describe sensitivity to these placements, not confidence limits.
#'
#' Shared cells do not establish paired samples, identical physical support,
#' comparable geology/era/reporting design, analytical agreement or independence.
#' Coordinates may be rounded/reused. This audit fits no concentration model,
#' estimates no method bias, supplies no p-value, and does not correct preferential
#' sampling or describe unsampled terrain. Spatial confidence limits remain
#' experimental. Review overlap before interpreting method-stratified estimates.
#' @seealso [analyze_elements()], [sampling_support_audit()], [declustering_sensitivity()]
#' @export
method_spatial_overlap <- function(df, cell_sizes, x_col, y_col,
                                   n_offsets = 2, method_col = "amethod") {
  annotated <- flag_analytical_methods(df, method_col)
  coordinates <- c(x_col, y_col)
  if (!is.character(x_col) || length(x_col) != 1L || is.na(x_col) || !nzchar(x_col) ||
      !is.character(y_col) || length(y_col) != 1L || is.na(y_col) || !nzchar(y_col) ||
      x_col == y_col || !all(coordinates %in% names(df)) ||
      !is.numeric(df[[x_col]]) || !is.numeric(df[[y_col]])) {
    stop("Supply distinct numeric projected coordinate columns")
  }
  if (method_col %in% coordinates) stop("Coordinate columns must differ from the method column")
  if (!is.numeric(cell_sizes) || !length(cell_sizes) ||
      any(!is.finite(cell_sizes)) || any(cell_sizes <= 0)) {
    stop("cell_sizes must contain positive finite widths")
  }
  cell_sizes <- sort(unique(cell_sizes))
  if (!is.numeric(n_offsets) || length(n_offsets) != 1L || !is.finite(n_offsets) ||
      n_offsets < 1 || n_offsets != floor(n_offsets) || n_offsets^2 > .Machine$integer.max) {
    stop("n_offsets must be a positive whole number with representable squared count")
  }
  if ("species" %in% names(df)) {
    species <- trimws(as.character(df$species))
    species <- unique(species[!is.na(species) & nzchar(species)])
    if (length(species) > 1L) stop("Select one element/species before auditing method overlap")
  }
  intervals <- .censored_intervals(df, "detected_val", "left_cen", "left_limit",
    "right_cen", "right_limit", "method_spatial_overlap")
  spatial <- intervals$ok & is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  anchor_x <- if (any(spatial)) min(df[[x_col]][spatial]) else NA_real_
  anchor_y <- if (any(spatial)) min(df[[y_col]][spatial]) else NA_real_
  codes <- .analytical_method_codes(df, method_col)
  levels <- sort(unique(codes), na.last = TRUE)
  index <- match(codes, levels)
  method_rows <- lapply(seq_along(levels), function(i) {
    rows <- index == i
    use <- rows & spatial
    first <- which(rows)[1L]
    tibble::tibble(method_code = levels[i],
      method_information = annotated$method_review_status[first],
      method_review_source = annotated$method_review_source[first],
      n_input = sum(rows), n_usable_chemistry = sum(rows & intervals$ok),
      n_excluded_chemistry = sum(rows & !intervals$ok),
      n_spatial = sum(use), n_unlocated_usable = sum(rows & intervals$ok & !spatial),
      n_coordinate_pairs = nrow(unique(data.frame(x = df[[x_col]][use], y = df[[y_col]][use]))))
  })
  method_counts <- if (length(method_rows)) dplyr::bind_rows(method_rows) else tibble::tibble(
    method_code = character(), method_information = character(), method_review_source = character(),
    n_input = integer(), n_usable_chemistry = integer(), n_excluded_chemistry = integer(),
    n_spatial = integer(), n_unlocated_usable = integer(), n_coordinate_pairs = integer())
  known <- which(!is.na(levels) & levels != "UNKNOWN")
  pairs <- if (length(known) >= 2L) utils::combn(known, 2) else matrix(integer(), nrow = 2)
  n_pairs <- ncol(pairs)
  n_origins <- n_offsets^2
  offsets <- (seq_len(n_offsets) - 1) / n_offsets
  summary_range <- function(values) {
    if (all(is.na(values))) c(NA_real_, NA_real_, NA_real_) else
      c(min(values, na.rm = TRUE), mean(values, na.rm = TRUE), max(values, na.rm = TRUE))
  }
  method_results <- pair_results <- vector("list", length(cell_sizes))
  for (s in seq_along(cell_sizes)) {
    size <- cell_sizes[s]
    x <- (df[[x_col]][spatial] - anchor_x) / size
    y <- (df[[y_col]][spatial] - anchor_y) / size
    if (any(!is.finite(x)) || any(!is.finite(y)) || any(x > 2^52) || any(y > 2^52)) {
      stop("cell_size is too small relative to the coordinate extent")
    }
    occupied <- matrix(NA_real_, nrow = length(levels), ncol = n_origins)
    shared <- union <- jaccard <- matrix(NA_real_, nrow = n_pairs, ncol = n_origins)
    origin <- 0L
    for (dx in offsets) for (dy in offsets) {
      origin <- origin + 1L
      ix <- floor(x + dx); iy <- floor(y + dy)
      # Use exact integer keys where possible; integer text must not round adjacent cells.
      width <- if (length(ix)) max(ix) + 1 else 1
      keys <- if (!length(iy) || width * (max(iy) + 1) <= 2^52) ix + width * iy else
        paste(sprintf("%.0f", ix), sprintf("%.0f", iy), sep = ":")
      cell_sets <- lapply(seq_along(levels), function(i) unique(keys[index[spatial] == i]))
      present <- which(method_counts$n_spatial > 0L)
      occupied[present, origin] <- lengths(cell_sets)[present]
      for (p in seq_len(n_pairs)) {
        a <- pairs[1L, p]; b <- pairs[2L, p]
        if (!length(cell_sets[[a]]) || !length(cell_sets[[b]])) next
        shared[p, origin] <- length(intersect(cell_sets[[a]], cell_sets[[b]]))
        union[p, origin] <- length(cell_sets[[a]]) + length(cell_sets[[b]]) - shared[p, origin]
        jaccard[p, origin] <- shared[p, origin] / union[p, origin]
      }
    }
    methods <- method_counts
    methods$cell_size <- rep(size, nrow(methods))
    methods$n_origins <- rep(n_origins, nrow(methods))
    methods$occupied_cells_min <- vapply(seq_along(levels), function(i) summary_range(occupied[i, ])[1L], 0.0)
    methods$occupied_cells_mean <- vapply(seq_along(levels), function(i) summary_range(occupied[i, ])[2L], 0.0)
    methods$occupied_cells_max <- vapply(seq_along(levels), function(i) summary_range(occupied[i, ])[3L], 0.0)
    method_results[[s]] <- methods
    pair_table <- tibble::tibble(cell_size = rep(size, n_pairs),
      method_a = levels[pairs[1L, ]], method_b = levels[pairs[2L, ]],
      n_spatial_a = method_counts$n_spatial[pairs[1L, ]], n_spatial_b = method_counts$n_spatial[pairs[2L, ]],
      n_origins = rep(n_origins, n_pairs),
      n_origins_evaluated = rowSums(!is.na(jaccard)),
      overlap_status = rep("unavailable_spatial_data", n_pairs))
    pair_table$overlap_status[pair_table$n_origins_evaluated > 0L] <- "descriptive_overlap_available"
    for (metric in c("shared_cells", "union_cells", "jaccard")) {
      values <- switch(metric, shared_cells = shared, union_cells = union, jaccard = jaccard)
      for (k in seq_along(c("min", "mean", "max"))) {
        pair_table[[paste0(metric, "_", c("min", "mean", "max")[k])]] <-
          vapply(seq_len(n_pairs), function(p) summary_range(values[p, ])[k], 0.0)
      }
    }
    pair_table$fraction_origins_with_overlap <- vapply(seq_len(n_pairs), function(p) {
      if (all(is.na(shared[p, ]))) NA_real_ else mean(shared[p, ] > 0, na.rm = TRUE)
    }, 0.0)
    pair_results[[s]] <- pair_table
  }
  list(support = tibble::tibble(n_input = nrow(df), n_usable_chemistry = sum(intervals$ok),
    n_excluded_chemistry = sum(!intervals$ok), n_spatial = sum(spatial),
    n_unlocated_usable = sum(intervals$ok & !spatial),
    n_named_codes_input = length(known), n_named_codes_spatial = sum(method_counts$n_spatial[known] > 0L),
    n_unresolved_spatial = sum(spatial & (is.na(codes) | codes == "UNKNOWN")),
    grid_anchor_x = anchor_x, grid_anchor_y = anchor_y),
    methods = dplyr::bind_rows(method_results), pairs = dplyr::bind_rows(pair_results))
}
