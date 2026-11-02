#' Flag analytical methods requiring measurement review
#'
#' Annotate reported observations without changing values, censoring fields,
#' validation status or row counts. This is a limited USGS method-code review,
#' not a certification of quantitative precision or statistical validity.
#'
#' @param df A data frame of reported observations.
#' @param method_col Method-code column. If absent, missing method information
#'   is reported rather than inferred from another field.
#' @return The input as a tibble with `method_review_status`,
#'   `review_semiquantitative` and `method_review_source`. Existing columns with
#'   those derived names are recalculated. Original method labels are preserved.
#'   The logical review flag is TRUE for documented semiquantitative codes and
#'   NA otherwise; NA does not mean that the method is quantitative.
#' @details
#' Currently recognizes USGS `ES_SQ` and `ES_H2O_P` as semiquantitative methods.
#' Matching ignores letter case and surrounding whitespace. Missing/blank codes
#' receive `method_not_supplied`, the literal code `UNKNOWN` receives
#' `method_unspecified`, and other codes receive `method_not_reviewed`.
#' Recognized codes receive `review_semiquantitative` with a documentation URL.
#' The not-reviewed status is no decision from this limited rule, not a finding
#' that another analytical method is faulty.
#'
#' A method family does not establish a unique reporting-class table or a
#' record-specific detection-limit design. No numeric class bounds, laboratory,
#' reporting era or missing method are invented. Unqualified values under these
#' methods need not represent exact continuous measurements; uncensored counts
#' elsewhere describe the input coding, not measurement precision. Estimates
#' are not filtered or changed, and clearing count-review flags does not resolve
#' method concerns. Spatial confidence limits remain experimental.
#'
#' @references
#' USGS Data Series 759, analytical-method code definitions:
#' \url{https://pubs.usgs.gov/ds/759/contents/DS759_pamphlet.pdf}.
#' @export
flag_analytical_methods <- function(df, method_col = "amethod") {
  if (!is.data.frame(df)) stop("df must be a data frame")
  reserved <- c("method_review_status", "review_semiquantitative", "method_review_source")
  if (!is.character(method_col) || length(method_col) != 1L || is.na(method_col) ||
      !nzchar(trimws(method_col)) || method_col %in% reserved) {
    stop("method_col must name one method column, not a derived review column")
  }
  code <- .analytical_method_codes(df, method_col)
  missing <- is.na(code)
  unspecified <- !missing & code == "UNKNOWN"
  semi <- !missing & code %in% c("ES_SQ", "ES_H2O_P")
  out <- tibble::as_tibble(df)
  out$method_review_status <- rep("method_not_reviewed", nrow(out))
  out$method_review_status[missing] <- "method_not_supplied"
  out$method_review_status[unspecified] <- "method_unspecified"
  out$method_review_status[semi] <- "review_semiquantitative"
  out$review_semiquantitative <- rep(NA, nrow(out))
  out$review_semiquantitative[semi] <- TRUE
  out$method_review_source <- rep(NA_character_, nrow(out))
  out$method_review_source[semi] <- "https://pubs.usgs.gov/ds/759/contents/DS759_pamphlet.pdf"
  out
}

.analytical_method_codes <- function(df, method_col = "amethod") {
  if (!method_col %in% names(df)) return(rep(NA_character_, nrow(df)))
  values <- df[[method_col]]
  if (!is.atomic(values) || !is.null(dim(values))) stop("Method codes must be an atomic vector")
  code <- toupper(trimws(as.character(values)))
  code[!is.na(code) & !nzchar(code)] <- NA_character_
  code
}

# Method counts use the same interval-validity mask as the estimator, with an
# explicit caller-supplied subset when coordinates or weights restrict a summary.
.summarize_analytical_methods <- function(df, usable = NULL, intervals = NULL, include_methods = TRUE) {
  annotated <- flag_analytical_methods(df)
  # Reuse the caller's validated intervals for custom column names and weights.
  if (is.null(intervals)) intervals <- .censored_intervals(df, "detected_val", "left_cen", "left_limit",
    "right_cen", "right_limit", "analytical method summary")
  coordinate_subset <- !is.null(usable)
  if (is.null(usable)) usable <- intervals$ok
  if (!is.logical(usable) || length(usable) != nrow(df) || anyNA(usable) ||
      any(usable & !intervals$ok)) stop("Invalid usable-observation mask for method summary")
  semi <- annotated$review_semiquantitative %in% TRUE
  used <- sum(usable)
  semi_used <- sum(semi & usable)
  summary <- tibble::tibble(
    n_semiquantitative_input = sum(semi), n_semiquantitative_used = semi_used,
    n_method_unspecified_used = sum(usable & annotated$method_review_status %in%
      c("method_not_supplied", "method_unspecified")),
    n_method_not_reviewed_used = sum(usable & annotated$method_review_status == "method_not_reviewed"),
    review_semiquantitative = if (semi_used > 0L) TRUE else NA,
    method_review_status = if (used == 0L) "no_usable_data" else if (semi_used > 0L)
      "review_semiquantitative" else "method_information_unresolved"
  )
  if (!include_methods) return(list(summary = summary))
  details <- tibble::tibble(
    method_code = .analytical_method_codes(df),
    method_review_status = annotated$method_review_status,
    review_semiquantitative = annotated$review_semiquantitative,
    method_review_source = annotated$method_review_source,
    usable = usable, uncensored = usable & intervals$lower == intervals$upper,
    left = usable & intervals$lower == 0,
    right = usable & is.infinite(intervals$upper)
  )
  table <- dplyr::summarise(dplyr::group_by(details, .data$method_code,
    .data$method_review_status, .data$review_semiquantitative, .data$method_review_source),
    n_input = dplyr::n(), n_usable = sum(.data$usable), n_uncensored = sum(.data$uncensored),
    n_left = sum(.data$left), n_right = sum(.data$right), .groups = "drop")
  table$analysis_scope <- rep(if (coordinate_subset) "chemistry_and_coordinates" else
    "unweighted_estimate", nrow(table))
  list(summary = summary, methods = table)
}

# Saved plots must retain unresolved measurement concerns. Older empirical-CDF
# objects predate method counts; their missing metadata cannot be reconstructed.
.method_review_caption <- function(summary, n_total) {
  required <- c("n_semiquantitative_used", "n_method_unspecified_used", "n_method_not_reviewed_used")
  if (!all(required %in% names(summary))) {
    return("Analytical-method review unavailable in this saved result; measurement precision unresolved.")
  }
  paste(
    sprintf("Known semiquantitative codes: %d/%d usable records; measurement precision unresolved.",
      summary$n_semiquantitative_used, n_total),
    sprintf("Missing/unspecified methods: %d; other methods unreviewed: %d.",
      summary$n_method_unspecified_used, summary$n_method_not_reviewed_used),
    sep = "\n"
  )
}
