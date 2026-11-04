#' Audit recorded sampling support before pooling geochemical samples
#'
#' Describes missingness and differences in depth, horizon, size fraction,
#' sampling medium, or other supplied support fields. Optional exact-coordinate
#' summaries flag co-located records with different recorded support.
#'
#' @param df One metadata row per sample, or one element's sample table.
#' @param group_cols Optional analysis-domain columns, such as state.
#' @param support_cols Fields to inspect. NULL selects available columns from
#'   depth, horizon, mesh_size, primary_class, sample_source, method_collected,
#'   collection_protocol, and amethod. At least one field must be available.
#' @param id_col Sample identifier column; IDs must be nonblank and unique.
#' @param x_col,y_col Optional numeric coordinate column names. Both are needed
#'   for location summaries; finite coordinate pairs are compared exactly.
#' @param include_memberships If TRUE, require coordinates and append a
#'   per-record `memberships` table and global `coordinate_groups` table.
#'   Default FALSE preserves the four existing output components.
#' @return A list with fields (counts and missingness by field/domain),
#'   combinations (recorded support tuples and counts), locations (optional
#'   exact-coordinate counts and review flags), and available_fields (names used).
#'   Optional memberships preserve original columns, row order and raw IDs, adding
#'   recorded_position_group, n_records_at_position, position_membership_status and
#'   review_shared_position. Unlocated rows have missing group/count/review fields.
#'   coordinate_groups reports exact pairs across all analysis domains, support
#'   differences/gaps and n_analysis_domains. A pooled audit counts as one domain.
#' @details
#' This is a descriptive audit, not a determination that samples are comparable.
#' Blank strings are treated as missing; otherwise recorded labels are retained
#' after trimming whitespace. Differences in spelling or units are not reconciled.
#' Mixed recorded support means at least one inspected field has more than one
#' nonmissing value at a location. Missing fields independently trigger a support
#' gap flag; missing information does not establish equivalence.
#'
#' Exact coordinate matches may reflect rounded or reused coordinates, not the
#' same physical sampling support. This function does not delete samples,
#' deduplicate co-located records, infer spatial independence, or assign weights.
#' It is not a near-neighbour clustering method. Inspect coordinate quality and
#' analytical methods alongside these summaries before defining pooled domains.
#'
#' Optional memberships use exact numeric pairs without rounding or text keys.
#' Groups span analysis domains so a shared position is not silently split by a
#' state or other label. The original locations component remains domain-specific.
#' Position labels are deterministic within a cohort under row reordering, but
#' are not persistent site identifiers across cohort changes. Coordinates must
#' share a suitable reference system; none is inferred or transformed here.
#' Missing/nonfinite coordinates never form a shared unknown-position group.
#' review_shared_position is TRUE for repeated positions and NA otherwise;
#' singletons are not certified independent. Records are neither collapsed nor
#' assigned weights. Original columns must not use optional derived field names.
#' @export
sampling_support_audit <- function(df, group_cols = character(), support_cols = NULL,
                                    id_col = "lab_id", x_col = NULL, y_col = NULL,
                                    include_memberships = FALSE) {
  if (!is.logical(include_memberships) || length(include_memberships) != 1L || is.na(include_memberships)) {
    stop("include_memberships must be TRUE or FALSE")
  }
  derived <- c("recorded_position_group", "n_records_at_position", "position_membership_status",
               "review_shared_position", "n_analysis_domains")
  if (include_memberships && any(derived %in% names(df))) stop("Optional membership fields already exist; rename original columns")
  if (include_memberships && (is.null(x_col) || is.null(y_col))) stop("include_memberships requires both coordinate columns")
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) || !id_col %in% names(df)) {
    stop("Missing sample identifier column")
  }
  defaults <- c("depth", "horizon", "mesh_size", "primary_class", "sample_source",
                "method_collected", "collection_protocol", "amethod")
  if (is.null(support_cols)) support_cols <- intersect(defaults, names(df))
  valid_names <- function(x) is.character(x) && !anyNA(x) && !anyDuplicated(x) && all(nzchar(x))
  if (!valid_names(group_cols) || !valid_names(support_cols) || !length(support_cols) ||
      !all(c(group_cols, support_cols) %in% names(df))) stop("Missing or invalid audit columns")
  reserved <- c("field", "n_samples", "n_recorded", "n_missing", "n_values",
                "fraction_recorded", "mixed_recorded_support", "has_support_gaps")
  if (length(intersect(group_cols, support_cols)) || any(c(group_cols, support_cols) %in% reserved) ||
      id_col %in% c(group_cols, support_cols)) stop("Audit columns must be distinct from identifiers and output fields")
  ids <- trimws(as.character(df[[id_col]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Sample identifiers must be nonblank and unique")
  data <- tibble::as_tibble(df)
  for (column in support_cols) {
    value <- trimws(as.character(data[[column]]))
    value[!is.na(value) & !nzchar(value)] <- NA_character_
    data[[column]] <- value
  }
  grouped <- dplyr::group_by(data, dplyr::across(dplyr::all_of(group_cols)))
  fields <- dplyr::bind_rows(lapply(support_cols, function(column) {
    out <- dplyr::summarise(grouped, n_samples = dplyr::n(),
      n_recorded = sum(!is.na(.data[[column]])), n_missing = sum(is.na(.data[[column]])),
      n_values = dplyr::n_distinct(.data[[column]], na.rm = TRUE), .groups = "drop")
    out$field <- column
    out$fraction_recorded <- ifelse(out$n_samples > 0, out$n_recorded / out$n_samples, NA_real_)
    out
  }))
  combinations <- dplyr::count(data, dplyr::across(dplyr::all_of(c(group_cols, support_cols))), name = "n_samples")
  locations <- tibble::tibble()
  if (xor(is.null(x_col), is.null(y_col))) stop("Supply both coordinate columns")
  if (!is.null(x_col)) {
    if (!valid_names(c(x_col, y_col)) || length(x_col) != 1L || length(y_col) != 1L ||
        !all(c(x_col, y_col) %in% names(data)) ||
        !is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]]) ||
        is.complex(data[[x_col]]) || is.complex(data[[y_col]]) ||
        !is.null(dim(data[[x_col]])) || !is.null(dim(data[[y_col]]))) stop("Coordinates must be distinct numeric columns")
    if (any(c(x_col, y_col) %in% c(group_cols, support_cols, reserved, id_col))) stop("Coordinate columns overlap audit fields")
    finite <- is.finite(data[[x_col]]) & is.finite(data[[y_col]])
    locations <- data[finite, , drop = FALSE] %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, x_col, y_col)))) %>%
      dplyr::summarise(n_samples = dplyr::n(),
        mixed_recorded_support = any(vapply(dplyr::pick(dplyr::all_of(support_cols)),
          function(z) dplyr::n_distinct(z, na.rm = TRUE) > 1L, logical(1))),
        has_support_gaps = any(is.na(as.matrix(dplyr::pick(dplyr::all_of(support_cols))))),
        .groups = "drop")
  }
  out <- list(fields = fields, combinations = combinations, locations = locations,
              available_fields = support_cols)
  if (include_memberships) {
    pairs <- tibble::tibble(.x = data[[x_col]], .y = data[[y_col]])
    unique_pairs <- dplyr::arrange(dplyr::distinct(pairs[finite, , drop = FALSE]), .data$.x, .data$.y)
    unique_pairs$.index <- seq_len(nrow(unique_pairs))
    index <- dplyr::left_join(pairs, unique_pairs, by = c(".x", ".y"))$.index
    index[!finite] <- NA_integer_
    labels <- sprintf("position_%06d", seq_len(nrow(unique_pairs)))
    counts <- tabulate(index[finite], nbins = nrow(unique_pairs))
    memberships <- tibble::as_tibble(df)
    memberships$recorded_position_group <- labels[index]
    memberships$n_records_at_position <- counts[index]
    memberships$position_membership_status <- ifelse(finite, "recorded_position_group", "unavailable_coordinates")
    memberships$review_shared_position <- rep(NA, nrow(data))
    shared <- which(!is.na(index) & memberships$n_records_at_position > 1L)
    memberships$review_shared_position[shared] <- TRUE

    # Choose a temporary field name without overwriting any supplied column.
    temporary <- "..geocensor_position_index"
    while (temporary %in% names(data)) temporary <- paste0(temporary, "_")
    indexed <- data[finite, unique(c(x_col, y_col, group_cols, support_cols)), drop = FALSE]
    indexed[[temporary]] <- index[finite]
    global <- dplyr::summarise(
      dplyr::group_by(indexed, dplyr::across(dplyr::all_of(c(temporary, x_col, y_col)))),
      n_records_at_position = dplyr::n(),
      mixed_recorded_support = any(vapply(dplyr::pick(dplyr::all_of(support_cols)),
        function(z) dplyr::n_distinct(z, na.rm = TRUE) > 1L, logical(1))),
      has_support_gaps = any(is.na(as.matrix(dplyr::pick(dplyr::all_of(support_cols))))),
      n_analysis_domains = if (length(group_cols)) nrow(dplyr::distinct(dplyr::pick(dplyr::all_of(group_cols)))) else 1L,
      .groups = "drop")
    global$recorded_position_group <- labels[global[[temporary]]]
    global[[temporary]] <- NULL
    global$review_shared_position <- rep(NA, nrow(global))
    global$review_shared_position[global$n_records_at_position > 1L] <- TRUE
    out$memberships <- memberships
    out$coordinate_groups <- global[c("recorded_position_group", x_col, y_col,
      "n_records_at_position", "mixed_recorded_support", "has_support_gaps",
      "n_analysis_domains", "review_shared_position")]
  }
  out
}
