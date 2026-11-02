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
#' @return A list with fields (counts and missingness by field/domain),
#'   combinations (recorded support tuples and counts), locations (optional
#'   exact-coordinate counts and review flags), and available_fields (names used).
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
#' @export
sampling_support_audit <- function(df, group_cols = character(), support_cols = NULL,
                                    id_col = "lab_id", x_col = NULL, y_col = NULL) {
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
        !is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]])) stop("Coordinates must be distinct numeric columns")
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
  list(fields = fields, combinations = combinations, locations = locations,
       available_fields = support_cols)
}
