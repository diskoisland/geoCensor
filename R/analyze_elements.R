#' Analyze multiple USGS geochemical elements
#'
#' Runs validation, censored lognormal estimation, distribution comparisons,
#' and optional cell-declustering sensitivity separately for each element and
#' analysis group. All concentrations are converted to ug/kg.
#'
#' @param chem_df Long chemistry table with lab_id, species, qvalue, and unit.
#' @param elements Element/species names to analyze. NULL selects all nonblank
#'   species in the input. Requested species with no rows are reported as
#'   insufficient data when no grouping columns are specified.
#' @param main_df Optional sample metadata table with unique nonblank lab_id.
#' @param group_cols Optional column names defining separate analysis domains,
#'   for example state or sampling medium. Columns may be in chemistry or metadata.
#'   Missing group values are retained as a separate group. An empty vector
#'   analyzes each element over the entire supplied dataset.
#' @param cell_sizes Optional vector of projected cell widths. NULL skips
#'   declustering. Weights are recalculated for each element/group's usable rows.
#' @param x_col,y_col Projected coordinate column names, required for declustering.
#' @param n_offsets Grid origins per axis for declustering.
#' @param min_exact_review,min_exact_fraction_review Configurable count-based
#'   review thresholds passed to [flag_estimate_information()]. These flags do
#'   not remove estimates or establish statistical reliability.
#' @return A list of tibbles: estimates (input/exclusion counts and fit results),
#'   distributions (unweighted model comparisons), validation (counts by reason),
#'   declustering (optional sensitivity results), and joins (optional join audits).
#' @details
#' Duplicate sample IDs within an element and missing sample IDs raise an error;
#' select one result per sample and element before analysis. Unsupported units
#' and invalid chemistry are excluded and counted, rather than stopping other
#' elements. A species with no usable observations still appears in estimates.
#' No automatic element selection by significance or fitted concentration occurs.
#' Supply one medium/domain at a time or use group_cols to avoid pooling unlike
#' sample populations. Weighted confidence limits are not computed.
#' @export
analyze_elements <- function(chem_df, elements = NULL, main_df = NULL,
                             group_cols = character(), cell_sizes = NULL,
                             x_col = NULL, y_col = NULL, n_offsets = 10,
                             min_exact_review = 10, min_exact_fraction_review = 0.05) {
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(chem_df))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(chem_df)), collapse = ", "))
  }
  if (is.null(elements)) {
    elements <- sort(unique(as.character(chem_df$species)))
    elements <- elements[!is.na(elements) & nzchar(trimws(elements))]
  }
  if (!is.character(elements) || !length(elements) || anyNA(elements) ||
      any(!nzchar(trimws(elements))) || anyDuplicated(elements)) {
    stop("elements must contain unique nonblank species names")
  }
  if (!is.character(group_cols) || anyNA(group_cols) || anyDuplicated(group_cols) ||
      any(!nzchar(group_cols)) || any(group_cols %in% c("element", "species", "lab_id"))) {
    stop("group_cols must contain unique analysis-domain column names")
  }
  flag_estimate_information(tibble::tibble(n_total = integer(), n_left = integer(),
    n_right = integer(), geomean = double()), min_exact_review, min_exact_fraction_review)
  spatial <- !is.null(cell_sizes)
  if (spatial && (is.null(x_col) || is.null(y_col))) {
    stop("x_col and y_col are required for declustering; use projected coordinates")
  }
  # Split once so repeated analyses do not repeatedly scan the full chemistry table.
  pieces <- split(chem_df, factor(as.character(chem_df$species), levels = elements))
  results <- lapply(elements, function(element) {
    observations <- build_element(pieces[[element]], element)
    ids <- trimws(as.character(observations$lab_id))
    if (anyNA(ids) || any(ids == "")) stop(element, ": missing or blank lab_id")
    if (anyDuplicated(ids)) stop(element, ": duplicate lab_id; select one result per sample and element")
    observations$lab_id <- ids
    audit <- tibble::tibble()
    if (!is.null(main_df)) {
      extras <- setdiff(unique(c(group_cols, x_col, y_col)), names(observations))
      observations <- join_soil_main(observations, main_df, extra_cols = extras)
      audit <- dplyr::bind_cols(tibble::tibble(element = element), attr(observations, "join_audit"))
    }
    if (!all(group_cols %in% names(observations))) stop("Missing grouping columns")
    if (spatial) {
      # Validate spatial arguments before splitting into possibly empty groups.
      declustering_sensitivity(observations[FALSE, ], cell_sizes, x_col, y_col, n_offsets)
    }
    groups <- if (length(group_cols)) {
      dplyr::group_split(dplyr::group_by(observations, dplyr::across(dplyr::all_of(group_cols))))
    } else list(observations)
    analyzed <- lapply(groups, function(x) {
      key <- tibble::tibble(element = element)
      if (length(group_cols)) key <- dplyr::bind_cols(key, x[1, group_cols, drop = FALSE])
      add_key <- function(tbl) dplyr::bind_cols(key[rep(1L, nrow(tbl)), , drop = FALSE], tbl)
      fit <- flag_estimate_information(fit_cen_lognorm(x),
        min_exact_review, min_exact_fraction_review)
      counts <- tibble::tibble(n_input = nrow(x), n_excluded = nrow(x) - fit$n_total,
                               value_unit = "ug/kg")
      list(
        estimates = add_key(dplyr::bind_cols(counts, fit)),
        distributions = add_key(compare_censored_distributions(x)),
        validation = add_key(dplyr::count(x, .data$validation_status, name = "n_observations")),
        declustering = if (spatial) add_key(declustering_sensitivity(x, cell_sizes, x_col, y_col, n_offsets)) else tibble::tibble()
      )
    })
    out <- lapply(c("estimates", "distributions", "validation", "declustering"), function(name) {
      dplyr::bind_rows(lapply(analyzed, `[[`, name))
    })
    names(out) <- c("estimates", "distributions", "validation", "declustering")
    out$joins <- audit
    out
  })
  names_out <- c("estimates", "distributions", "validation", "declustering", "joins")
  stats::setNames(lapply(names_out, function(name) dplyr::bind_rows(lapply(results, `[[`, name))), names_out)
}
