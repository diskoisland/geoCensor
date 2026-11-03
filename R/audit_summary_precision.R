#' Audit source-summary precision without replacing reported values
#'
#' Compare a reported value with the first value/method entry in its USGS-style
#' summary. Source-value candidates require explicit matching method codes and
#' compatible usable censoring encodings. No input value or row is changed.
#'
#' @param df A data frame containing reported chemistry.
#' @param value_col Name of the reported value column; it must exist.
#' @param method_col Name of the method-code column. Missing codes are reported.
#' @param summary_col Name of the semicolon-separated value/method summary column.
#' @param qualifier_col Optional qualifier column name; absent columns are treated
#'   as unspecified qualifiers by [parse_usgs_q()].
#' @return The input as a tibble with eight audit columns:
#'   `summary_first_value_raw`, `summary_first_method_raw`,
#'   `summary_method_matches`, `summary_reported_coding_status`,
#'   `summary_entry_coding_status`, `summary_qvalue_candidate`,
#'   `summary_numeric_delta` and `summary_precision_status`.
#' @details
#' Only the first summary entry is inspected; it must contain one numeric value
#' and one method separated by a comma. Method matching ignores case and surrounding
#' whitespace. Blank, missing and literal `UNKNOWN` method labels do not establish
#' method association, even if repeated in both fields. Missing method information
#' is never filled from the summary.
#'
#' A candidate is review evidence, not an automatically recovered measurement or
#' cleared model input. `different_reported_value` and `different_reported_limit`
#' expose source discrepancies without selecting either representation. An explicit
#' upper-censoring suffix that rounds to the reported value receives
#' `summary_upper_limit_evidence`. Other censoring-class conflicts expose no
#' candidate. Unknown qualifiers and nonpositive uncensored reported values remain
#' unusable according to [parse_usgs_q()].
#'
#' `summary_numeric_delta` is the summary concentration or censoring limit minus
#' its reported counterpart, in unchanged input units. For upper-limit evidence
#' on an otherwise unqualified reported value, the numeric reported value is used
#' as the counterpart. The candidate retains the summary's encoded numeric value,
#' including a negative lower-limit code or explicit upper suffix.
#'
#' This audit does not establish analytical accuracy, measurement precision,
#' specimen/site identity, digest recovery, sampling independence or confidence
#' coverage. It does not alter existing parser or estimator behavior. Derived
#' audit columns must not already exist in the input. Spatial limits remain
#' experimental.
#' @export
audit_summary_precision <- function(df, value_col = "qvalue", method_col = "amethod",
                                    summary_col = "summary", qualifier_col = "qualifier") {
  if (!is.data.frame(df)) stop("df must be a data frame")
  columns <- list(value_col, method_col, summary_col, qualifier_col)
  valid <- vapply(columns, function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(trimws(x)), logical(1))
  if (!all(valid) || anyDuplicated(unlist(columns))) stop("Audit column names must be distinct nonempty strings")
  if (!value_col %in% names(df)) stop("Missing reported value column")
  reserved <- .summary_precision_columns()
  if (any(reserved %in% names(df))) stop("Derived summary audit columns already exist")
  for (column in intersect(unlist(columns), names(df))) {
    value <- df[[column]]
    if (!is.atomic(value) || !is.null(dim(value))) stop("Audit columns must be atomic vectors")
  }
  n <- nrow(df)
  get_text <- function(column) {
    if (!column %in% names(df)) return(rep(NA_character_, n))
    value <- trimws(as.character(df[[column]]))
    value[!is.na(value) & !nzchar(value)] <- NA_character_
    value
  }
  reported <- get_text(value_col)
  method <- toupper(get_text(method_col))
  summaries <- get_text(summary_col)
  qualifier <- get_text(qualifier_col)
  first_value <- first_method <- rep(NA_character_, n)
  entry_present <- !is.na(summaries)
  parsed_entry <- rep(FALSE, n)
  for (i in which(entry_present)) {
    entry <- trimws(sub(";.*$", "", summaries[i]))
    if (nchar(entry) - nchar(gsub(",", "", entry, fixed = TRUE)) != 1L) next
    tokens <- strsplit(entry, ",", fixed = TRUE)[[1L]]
    if (length(tokens) == 2L && all(nzchar(trimws(tokens)))) {
      first_value[i] <- trimws(tokens[1L])
      first_method[i] <- trimws(tokens[2L])
      parsed_entry[i] <- TRUE
    }
  }
  matches <- rep(NA, n)
  known_method <- !is.na(method) & method != "UNKNOWN" & !is.na(first_method) &
    toupper(first_method) != "UNKNOWN"
  matches[known_method] <- method[known_method] == toupper(first_method[known_method])
  raw <- parse_usgs_q(reported, qualifier)
  source <- parse_usgs_q(first_value)
  candidate <- delta <- rep(NA_real_, n)
  status <- rep("not_available_summary", n)
  status[entry_present] <- "unparsed_summary"
  status[parsed_entry] <- "method_information_missing"
  status[parsed_entry & !is.na(matches) & !matches] <- "method_mismatch"
  eligible <- parsed_entry & !is.na(matches) & matches
  status[eligible] <- "reported_coding_needs_review"
  raw_ok <- eligible & raw$validation_status == "usable"
  status[raw_ok] <- "summary_coding_needs_review"
  usable <- raw_ok & source$validation_status == "usable"
  raw_number <- suppressWarnings(as.numeric(reported))
  source_number <- suppressWarnings(as.numeric(first_value))
  for (i in which(usable)) {
    upper_evidence <- source$right_cen[i] && !raw$left_cen[i] &&
      grepl("\\.1{4,5}$", first_value[i]) &&
      abs(round(source_number[i], 2) - raw_number[i]) < 1e-8 &&
      (!raw$right_cen[i] || raw$right_limit[i] != source$right_limit[i])
    if (upper_evidence) {
      counterpart <- if (raw$right_cen[i]) raw$right_limit[i] else raw_number[i]
      candidate[i] <- source_number[i]
      delta[i] <- source$right_limit[i] - counterpart
      status[i] <- "summary_upper_limit_evidence"
    } else if (raw$left_cen[i] != source$left_cen[i] || raw$right_cen[i] != source$right_cen[i]) {
      status[i] <- "censoring_encoding_difference"
    } else {
      candidate[i] <- source_number[i]
      if (raw$left_cen[i]) {
        delta[i] <- source$left_limit[i] - raw$left_limit[i]
      } else if (raw$right_cen[i]) {
        delta[i] <- source$right_limit[i] - raw$right_limit[i]
      } else {
        delta[i] <- source$detected_val[i] - raw$detected_val[i]
      }
      status[i] <- if (delta[i] == 0) "same_reported_value" else if (raw$left_cen[i] || raw$right_cen[i])
        "different_reported_limit" else "different_reported_value"
    }
  }
  out <- tibble::as_tibble(df)
  out$summary_first_value_raw <- first_value
  out$summary_first_method_raw <- first_method
  out$summary_method_matches <- matches
  out$summary_reported_coding_status <- raw$validation_status
  out$summary_entry_coding_status <- ifelse(parsed_entry, source$validation_status, NA_character_)
  out$summary_qvalue_candidate <- candidate
  out$summary_numeric_delta <- delta
  out$summary_precision_status <- status
  out
}

.summary_precision_columns <- function() c(
  "summary_first_value_raw", "summary_first_method_raw", "summary_method_matches",
  "summary_reported_coding_status", "summary_entry_coding_status", "summary_qvalue_candidate",
  "summary_numeric_delta", "summary_precision_status")
