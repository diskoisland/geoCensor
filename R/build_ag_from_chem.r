#' Build Ag (or other element) censoring fields from a chemistry data frame
#'
#' @param chem_df Data frame with at least columns: lab_id, species, qvalue, and
#'   unit. A qualifier column is optional and treated as missing when absent.
#' @param element Element symbol to filter (e.g. "Ag").
#' @param to_ugkg Logical; if TRUE, convert ppm (mg/kg) → µg/kg.
#'
#' @return A data frame retaining each input row for that element, including
#'         detected_val, left_cen, left_limit, right_cen, right_limit, value_unit.
#' `validation_status` identifies unusable values and unsupported units.
#' Supported mass-fraction units are ppm, mg/kg, ppb, ug/kg (including micro
#' symbols), pct, and percent. Duplicate sample rows are retained.
#' If `summary` contains an explicit upper-censoring suffix in its first
#' (best-ranked) analysis, that value rounds to the reported qvalue, and a known
#' method matches with compatible usable coding, its censoring is recovered.
#' Missing/UNKNOWN methods, mismatches, malformed entries and coding conflicts
#' block recovery. The original qvalue/qualifier are then parsed as reported.
#' `summary_recovery_status` records recovery, unnecessary recovery, no upper
#' evidence, or a specific blocking reason; it does not certify measurement QA.
#' `qvalue_raw` retains the original value;
#' `censoring_source` records whether the limit came from summary or qvalue.
#' [flag_analytical_methods()] adds separate method-review fields. Method
#' concerns do not change `validation_status`, censoring fields or values;
#' semiquantitative reporting-class bounds are not inferred.
#' @export
build_ag_from_chem <- function(chem_df,
                               element = "Ag",
                               to_ugkg = TRUE) {
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(chem_df))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(chem_df)), collapse = ", "))
  }
  if (!is.character(element) || length(element) != 1L || is.na(element) || !nzchar(element)) {
    stop("element must be one nonempty element symbol")
  }
  if (!is.logical(to_ugkg) || length(to_ugkg) != 1L || is.na(to_ugkg)) {
    stop("to_ugkg must be TRUE or FALSE")
  }
  if (!"qualifier" %in% names(chem_df)) {
    chem_df <- dplyr::mutate(chem_df, qualifier = NA_character_)
  }

  df <- chem_df %>%
    dplyr::filter(.data$species == element) %>%
    dplyr::mutate(
      qvalue_raw = as.character(.data$qvalue),
      q_num      = suppressWarnings(as.numeric(.data$qvalue_raw)),
      unit_norm  = stringr::str_to_lower(stringr::str_trim(.data$unit))
    )
  
  effective <- .qvalue_for_censoring(df)
  cen <- parse_usgs_q(effective$value, df$qualifier)
  
  out <- dplyr::bind_cols(df, cen)
  out$censoring_source <- effective$source
  out$summary_recovery_status <- effective$status
  
  if (to_ugkg) {
    out <- out %>%
      dplyr::mutate(
        unit_factor  = .concentration_factor_ugkg(.data$unit),
        validation_status = dplyr::if_else(
          .data$validation_status == "usable" & !is.finite(.data$unit_factor),
          "unsupported_unit", .data$validation_status
        ),
        detected_val = .data$detected_val * .data$unit_factor,
        left_limit   = .data$left_limit * .data$unit_factor,
        right_limit  = .data$right_limit * .data$unit_factor,
        value_unit   = dplyr::if_else(is.finite(.data$unit_factor), "ug/kg", .data$unit)
      )
  } else {
    out <- out %>% dplyr::mutate(value_unit = .data$unit)
  }
  
  flag_analytical_methods(out)
}

.concentration_factor_ugkg <- function(unit) {
  normalized <- stringr::str_to_lower(stringr::str_trim(as.character(unit)))
  dplyr::case_when(
    normalized %in% c("ppm", "mg/kg") ~ 1000,
    normalized %in% c("ppb", "ug/kg", "\u00b5g/kg", "\u03bcg/kg") ~ 1,
    normalized %in% c("pct", "%", "percent") ~ 1e7,
    TRUE ~ NA_real_
  )
}

# Best-value CSV exports may round qvalue to two decimal places. The first
# summary value is the best-ranked original analysis (USGS data dictionary).
# Inspect only potential upper-suffix entries. The standalone audit supplies the
# same strict first-entry/method/coding rules without auditing millions of ordinary
# rows during routine parsing. Blocked evidence leaves original coding intact.
.qvalue_for_censoring <- function(df) {
  value <- as.character(df$qvalue)
  source <- rep("qvalue", nrow(df))
  status <- rep("no_summary", nrow(df))
  if ("summary" %in% names(df)) {
    summaries <- trimws(as.character(df$summary))
    supplied <- !is.na(summaries) & nzchar(summaries)
    status[supplied] <- "no_upper_limit_evidence"
    first <- trimws(sub(",.*$", "", sub(";.*$", "", summaries)))
    candidates <- which(supplied & grepl("\\.1{4,5}$", first))
    if (length(candidates)) {
      fields <- intersect(c("qvalue", "amethod", "summary", "qualifier"), names(df))
      checked <- audit_summary_precision(df[candidates, fields, drop = FALSE])
      labels <- c(
        summary_upper_limit_evidence = "recovered_summary_upper_limit",
        same_reported_value = "original_upper_limit_retained",
        different_reported_value = "blocked_value_disagreement",
        different_reported_limit = "blocked_limit_disagreement",
        censoring_encoding_difference = "blocked_censoring_conflict",
        method_mismatch = "blocked_method_mismatch",
        method_information_missing = "blocked_method_information_missing",
        reported_coding_needs_review = "blocked_reported_coding",
        summary_coding_needs_review = "blocked_summary_coding",
        unparsed_summary = "blocked_malformed_summary",
        not_available_summary = "no_summary")
      mapped <- unname(labels[checked$summary_precision_status])
      if (anyNA(mapped)) stop("Unrecognized summary-recovery review status")
      status[candidates] <- mapped
      recover <- checked$summary_precision_status == "summary_upper_limit_evidence"
      value[candidates[recover]] <- checked$summary_first_value_raw[recover]
      source[candidates[recover]] <- "summary_upper_limit"
    }
  }
  list(value = value, source = source, status = status)
}
