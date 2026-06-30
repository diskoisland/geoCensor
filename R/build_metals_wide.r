#' Build a wide table of metals (ppm → mg/kg) using substitution for censoring
#'
#' @param chem_df Chemistry long table (lab_id, species, qvalue, unit, qualifier).
#' @param elements Character vector of element symbols to include.
#' @param lod_factor Factor for left-censored substitution (default LOD/2).
#' @param suffix Suffix for column names (default "_mgkg").
#'
#' @return Wide data frame: one row per lab_id, one column per element.
#' @export
build_metals_wide <- function(chem_df,
                              elements,
                              lod_factor = 0.5,
                              suffix     = "_mgkg") {
  df <- chem_df %>%
    dplyr::filter(.data$species %in% elements) %>%
    dplyr::mutate(
      qvalue_raw = as.character(.data$qvalue),
      q_num      = readr::parse_number(.data$qvalue_raw),
      has_upper  = stringr::str_detect(.data$qvalue_raw, "\\.1111$")
    )
  
  left_cen   <- is.finite(df$q_num) & df$q_num < 0
  right_cen  <- df$has_upper
  left_limit <- ifelse(left_cen, abs(df$q_num), NA_real_)
  right_limit <- ifelse(
    right_cen,
    readr::parse_number(stringr::str_remove(df$qvalue_raw, "\\.1111$")),
    NA_real_
  )
  
  detected_val <- ifelse(left_cen | right_cen, NA_real_, df$q_num)
  
  detected_val <- ifelse(df$unit == "ppm", detected_val * 1000, detected_val)
  left_limit   <- ifelse(df$unit == "ppm", left_limit   * 1000, left_limit)
  right_limit  <- ifelse(df$unit == "ppm", right_limit  * 1000, right_limit)
  
  val_ugkg <- dplyr::case_when(
    left_cen  & is.finite(left_limit)  ~ lod_factor * left_limit,
    right_cen & is.finite(right_limit) ~ right_limit,
    TRUE                               ~ detected_val
  )
  val_mgkg <- val_ugkg / 1000
  
  df_sub <- df %>%
    dplyr::transmute(
      lab_id,
      species = paste0(.data$species, suffix),
      val_mgkg = val_mgkg
    )
  
  df_sub %>%
    tidyr::pivot_wider(
      id_cols     = .data$lab_id,
      names_from  = .data$species,
      values_from = .data$val_mgkg,
      values_fn   = list(val_mgkg = ~ mean(.x, na.rm = TRUE))
    )
}

#' Convenience: join metals wide + main + state filter
#'
#' @param chem_df Chemistry table.
#' @param main_df Main metadata (from read_clean_main).
#' @param elements Metals to include.
#' @param state Two-letter USPS code or full name; if NULL, keep all.
#'
#' @return Wide metals table joined with state (and other metadata).
#' @export
build_ag_metals_state <- function(chem_df,
                                  main_df,
                                  elements,
                                  state = NULL) {
  wide <- build_metals_wide(chem_df, elements = elements)
  
  out <- wide %>%
    dplyr::left_join(
      main_df %>% dplyr::select(.data$lab_id, .data$state, .data$lat_wgs84, .data$long_wgs84),
      by = "lab_id"
    )
  
  if (!is.null(state)) {
    state_code <- norm_state(state)
    out <- out %>% dplyr::filter(.data$state == state_code)
  }
  
  out
}