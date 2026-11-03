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
#' @param n_offsets Grid origins per axis for declustering and optional method
#'   overlap. The batch default is 10 (100 placements), rather than the
#'   standalone overlap audit's default of 2.
#' @param min_exact_review,min_exact_fraction_review Configurable count-based
#'   review thresholds passed to [flag_estimate_information()]. These flags do
#'   not remove estimates or establish statistical reliability.
#' @param stratify_methods If TRUE, append a `method_estimates` table of
#'   unweighted lognormal summaries by normalized `amethod` code within each
#'   element/domain. The default FALSE leaves the six existing components intact.
#' @param audit_method_overlap If TRUE, append method_overlap_support,
#'   method_coverage and method_overlap tables from [method_spatial_overlap()]
#'   within each element/domain. Requires cell_sizes and projected x_col/y_col.
#'   Independent of stratify_methods; default FALSE leaves existing outputs intact.
#' @param audit_source_precision If TRUE, append `source_precision` with the
#'   requested original chemistry rows and [audit_summary_precision()] fields,
#'   plus `source_precision_summary` counts by element/domain/status. Default
#'   FALSE leaves existing components intact. No candidate replaces model inputs.
#' @return A list of tibbles: estimates (input/exclusion counts and fit results),
#'   distributions (unweighted model comparisons), validation (counts by reason),
#'   declustering (optional sensitivity results), joins (optional join audits),
#'   and methods (per-code input and unweighted fit-used counts with provenance).
#'   Estimates also include separate method-review counts/status. Semiquantitative
#'   concerns do not invalidate records or alter numerical estimates.
#'   With stratify_methods = TRUE, a seventh component, method_estimates, adds
#'   method_code, method_information, method_review_source, analysis_scope,
#'   input/exclusion counts, estimates and the same information-review fields.
#'   With audit_method_overlap = TRUE, three additional components record data
#'   availability, coverage by method/scale and named-method pair overlap. They
#'   use sorted unique cell sizes and the same n_offsets value as declustering.
#'   With audit_source_precision = TRUE, source_precision retains requested input
#'   row order, raw IDs/values/units and all original chemistry columns. Grouping
#'   columns available only in metadata are appended using the validated ID join.
#'   Its candidates and numeric deltas remain in each row's original units;
#'   estimator outputs continue to use ug/kg. source_precision_summary counts all
#'   requested input rows, including those excluded from numerical estimation.
#' @details
#' Duplicate sample IDs within an element and missing sample IDs raise an error;
#' select one result per sample and element before analysis. Unsupported units
#' and invalid chemistry are excluded and counted, rather than stopping other
#' elements. A species with no usable observations still appears in estimates.
#' No automatic element selection by significance or fitted concentration occurs.
#' Supply one medium/domain at a time or use group_cols to avoid pooling unlike
#' sample populations. Weighted confidence limits are not computed.
#' Method counts in estimates use valid chemistry, while the methods table also
#' retains excluded input records. Declustering method counts instead use the
#' chemistry/coordinate subset. Uncensored counts describe input coding,
#' not verified measurement precision. Grouping names must not conflict with
#' result fields; rename a grouping column if necessary, or group by raw amethod.
#' Missing or unreviewed method information
#' is not interpreted as quantitative; see [flag_analytical_methods()].
#'
#' Optional method stratification trims whitespace and ignores code letter case,
#' without modifying raw input labels. Blank/absent codes form an explicitly
#' unresolved NA bucket; the literal UNKNOWN remains a separate unspecified
#' bucket. Such buckets do not establish a shared analytical procedure. All
#' recorded buckets are retained, including those without usable chemistry;
#' absent species do not create phantom method rows. Method summaries use valid
#' chemistry, irrespective of coordinate availability, and do not recalculate
#' spatial weights. Ordinary limits retain their independent-observation model
#' assumptions; spatial confidence limits remain experimental.
#'
#' Differences between strata are descriptive sensitivity results, not estimates
#' of analytical bias. Method families can also differ in geology, location,
#' sample support, era and reporting design. No method is selected by its fitted
#' concentration or significance; no missing method, calibration or reporting
#' interval is inferred. Review comparability before interpreting or pooling.
#' Optional overlap tables use the chemistry/coordinate-eligible subset while
#' retaining original input and unlocated/excluded counts. Missing comparisons
#' stay NA. Shared cells do not establish paired samples or analytical agreement;
#' grid-placement ranges are not confidence limits. Neither diagnostic modifies
#' existing estimates, chooses a method or changes declustering weights.
#'
#' Source-precision review compares the original reported value with the first
#' explicitly associated summary entry before any existing parser recovery or unit
#' conversion. It does not change that recovery, select a representation, certify
#' analytical reliability or filter observations. Missing summaries and unresolved
#' methods remain explicit. Requested species with no input rows produce no
#' phantom audit records or status counts. Grouping names must not conflict with
#' the eight source-audit fields or their status-count field `n_observations`.
#' @export
analyze_elements <- function(chem_df, elements = NULL, main_df = NULL,
                             group_cols = character(), cell_sizes = NULL,
                             x_col = NULL, y_col = NULL, n_offsets = 10,
                             min_exact_review = 10, min_exact_fraction_review = 0.05,
                             stratify_methods = FALSE, audit_method_overlap = FALSE,
                             audit_source_precision = FALSE) {
  if (!is.logical(stratify_methods) || length(stratify_methods) != 1L || is.na(stratify_methods)) {
    stop("stratify_methods must be TRUE or FALSE")
  }
  if (!is.logical(audit_method_overlap) || length(audit_method_overlap) != 1L || is.na(audit_method_overlap)) {
    stop("audit_method_overlap must be TRUE or FALSE")
  }
  if (!is.logical(audit_source_precision) || length(audit_source_precision) != 1L || is.na(audit_source_precision)) {
    stop("audit_source_precision must be TRUE or FALSE")
  }
  if (audit_method_overlap && is.null(cell_sizes)) stop("audit_method_overlap requires cell_sizes and projected coordinates")
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
  source_precision <- NULL
  if (audit_source_precision) {
    conflicts <- intersect(group_cols, c(.summary_precision_columns(), "n_observations"))
    if (length(conflicts)) stop("Grouping columns conflict with source-precision fields: ",
      paste(conflicts, collapse = ", "), ". Rename these grouping columns before analysis.", call. = FALSE)
    requested <- as.character(chem_df$species) %in% elements
    source_precision <- audit_summary_precision(chem_df[requested, , drop = FALSE])
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
      if (audit_method_overlap) method_spatial_overlap(observations[FALSE, ], cell_sizes, x_col, y_col, n_offsets)
    }
    groups <- if (length(group_cols)) {
      dplyr::group_split(dplyr::group_by(observations, dplyr::across(dplyr::all_of(group_cols))))
    } else list(observations)
    analyzed <- lapply(groups, function(x) {
      key <- tibble::tibble(element = element)
      if (length(group_cols)) key <- dplyr::bind_cols(key, x[1, group_cols, drop = FALSE])
      add_key <- function(tbl) {
        overlap <- intersect(names(key), names(tbl))
        if (length(overlap)) stop("Grouping columns conflict with result fields: ",
          paste(overlap, collapse = ", "), ". Rename these grouping columns before analysis.", call. = FALSE)
        dplyr::bind_cols(key[rep(1L, nrow(tbl)), , drop = FALSE], tbl)
      }
      fit <- flag_estimate_information(fit_cen_lognorm(x),
        min_exact_review, min_exact_fraction_review)
      method_audit <- .summarize_analytical_methods(x)
      # Keep batch review fields after the established count-review columns.
      method_fields <- names(method_audit$summary)
      fit <- fit[c(setdiff(names(fit), method_fields), method_fields)]
      counts <- tibble::tibble(n_input = nrow(x), n_excluded = nrow(x) - fit$n_total,
                               value_unit = "ug/kg")
      result <- list(
        estimates = add_key(dplyr::bind_cols(counts, fit)),
        distributions = add_key(compare_censored_distributions(x)),
        validation = add_key(dplyr::count(x, .data$validation_status, name = "n_observations")),
        methods = add_key(method_audit$methods),
        declustering = if (spatial) add_key(declustering_sensitivity(x, cell_sizes, x_col, y_col, n_offsets)) else tibble::tibble()
      )
      if (stratify_methods) result$method_estimates <- add_key(.method_stratified_estimates(
        x, method_fields, min_exact_review, min_exact_fraction_review))
      if (audit_method_overlap) {
        geographic <- method_spatial_overlap(x, cell_sizes, x_col, y_col, n_offsets)
        result$method_overlap_support <- add_key(geographic$support)
        result$method_coverage <- add_key(geographic$methods)
        result$method_overlap <- add_key(geographic$pairs)
      }
      result
    })
    components <- c("estimates", "distributions", "validation", "declustering", "methods",
      if (stratify_methods) "method_estimates",
      if (audit_method_overlap) c("method_overlap_support", "method_coverage", "method_overlap"))
    out <- lapply(components, function(name) {
      dplyr::bind_rows(lapply(analyzed, `[[`, name))
    })
    names(out) <- components
    out$joins <- audit
    out
  })
  names_out <- c("estimates", "distributions", "validation", "declustering", "joins", "methods",
    if (stratify_methods) "method_estimates",
    if (audit_method_overlap) c("method_overlap_support", "method_coverage", "method_overlap"))
  out <- stats::setNames(lapply(names_out, function(name) dplyr::bind_rows(lapply(results, `[[`, name))), names_out)
  if (audit_source_precision) {
    # The ordinary workflow above validates keys and joins. Attach only metadata
    # grouping fields here, preserving original source IDs and row order.
    metadata_groups <- setdiff(group_cols, names(source_precision))
    if (length(metadata_groups)) {
      matched <- match(trimws(as.character(source_precision$lab_id)),
        trimws(as.character(main_df$lab_id)))
      for (column in metadata_groups) source_precision[[column]] <- main_df[[column]][matched]
    }
    summary_input <- dplyr::bind_cols(
      tibble::tibble(element = as.character(source_precision$species)),
      source_precision[group_cols],
      tibble::tibble(summary_precision_status = source_precision$summary_precision_status))
    keys <- c("element", group_cols, "summary_precision_status")
    out$source_precision <- source_precision
    out$source_precision_summary <- dplyr::summarise(
      dplyr::group_by(summary_input, dplyr::across(dplyr::all_of(keys))),
      n_observations = dplyr::n(), .groups = "drop")
  }
  out
}

# These unweighted summaries partition the existing element/domain input, not
# the coordinate-eligible subset used for declustering. Unknown labels stay unknown.
.method_stratified_estimates <- function(df, method_fields, min_exact, min_exact_fraction) {
  one <- function(x, code) {
    information <- flag_analytical_methods(data.frame(amethod = code))
    fit <- flag_estimate_information(fit_cen_lognorm(x), min_exact, min_exact_fraction)
    fit <- fit[c(setdiff(names(fit), method_fields), method_fields)]
    dplyr::bind_cols(tibble::tibble(
      method_code = code, method_information = information$method_review_status,
      method_review_source = information$method_review_source,
      analysis_scope = "unweighted_method_stratum",
      n_input = nrow(x), n_excluded = nrow(x) - fit$n_total, value_unit = "ug/kg"), fit)
  }
  if (!nrow(df)) return(one(df, NA_character_)[FALSE, ])
  codes <- .analytical_method_codes(df)
  levels <- sort(unique(codes), na.last = TRUE)
  indices <- split(seq_len(nrow(df)), match(codes, levels))
  dplyr::bind_rows(lapply(seq_along(levels), function(i) one(df[indices[[as.character(i)]], , drop = FALSE], levels[i])))
}
