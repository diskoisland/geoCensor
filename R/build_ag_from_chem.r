#' Build Ag (or other element) censoring fields from a chemistry data frame
#'
#' @param chem_df Data frame with at least columns: lab_id, species, qvalue, unit, qualifier.
#' @param element Element symbol to filter (e.g. "Ag").
#' @param to_ugkg Logical; if TRUE, convert ppm (mg/kg) → µg/kg.
#'
#' @return A data frame with one row per lab_id for that element, including
#'         detected_val, left_cen, left_limit, right_cen, right_limit, value_unit.
#' @export
build_ag_from_chem <- function(chem_df,
                               element = "Ag",
                               to_ugkg = TRUE) {
  df <- chem_df %>%
    dplyr::filter(.data$species == element) %>%
    dplyr::mutate(
      qvalue_raw = as.character(.data$qvalue),
      q_num      = readr::parse_number(.data$qvalue_raw),
      unit_norm  = stringr::str_to_lower(stringr::str_trim(.data$unit))
    )
  
  cen <- parse_usgs_q(df$qvalue_raw, df$qualifier)
  
  out <- dplyr::bind_cols(df, cen)
  
  if (to_ugkg) {
    out <- out %>%
      dplyr::mutate(
        unit_factor  = dplyr::case_when(
          .data$unit_norm == "ppm" ~ 1000,
          .data$unit_norm == "pct" ~ 1e7,
          TRUE                     ~ NA_real_
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
