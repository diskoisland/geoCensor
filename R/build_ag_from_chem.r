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
#' (best-ranked) analysis and that value rounds to the reported qvalue, its
#' censoring information is recovered. `qvalue_raw` retains the original value;
#' `censoring_source` records whether the limit came from summary or qvalue.
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
  
  out
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
# Recover only explicit upper-censoring suffixes that agree with rounded qvalue.
.qvalue_for_censoring <- function(df) {
  value <- as.character(df$qvalue)
  source <- rep("qvalue", nrow(df))
  if ("summary" %in% names(df)) {
    first <- trimws(sub(",.*$", "", as.character(df$summary)))
    original <- suppressWarnings(as.numeric(first))
    reported <- suppressWarnings(as.numeric(value))
    recover <- grepl("\\.1{4,5}$", first) & is.finite(original) &
      is.finite(reported) & abs(round(original, 2) - reported) < 1e-8
    value[recover] <- first[recover]
    source[recover] <- "summary_upper_limit"
  }
  list(value = value, source = source)
}
