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
      q_num      = readr::parse_number(.data$qvalue_raw)
    )
  
  cen <- parse_usgs_q(df$qvalue_raw, df$qualifier)
  
  out <- dplyr::bind_cols(df, cen)
  
  if (to_ugkg) {
    out <- out %>%
      dplyr::mutate(
        detected_val = dplyr::if_else(.data$unit == "ppm", .data$detected_val * 1000, .data$detected_val),
        left_limit   = dplyr::if_else(.data$unit == "ppm", .data$left_limit   * 1000, .data$left_limit),
        right_limit  = dplyr::if_else(.data$unit == "ppm", .data$right_limit  * 1000, .data$right_limit),
        value_unit   = dplyr::if_else(.data$unit == "ppm", "ug/kg", .data$unit)
      )
  } else {
    out <- out %>% dplyr::mutate(value_unit = .data$unit)
  }
  
  out
}