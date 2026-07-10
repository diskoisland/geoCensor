#' Build a wide table of metals (mg/kg) using substitution for censoring
#'
#' Censoring is parsed with [parse_usgs_q()], so left/right censoring encoded via
#' negative qvalues, the `.1111` upper-limit suffix, **or** the `L`/`<`/`G`/`>`
#' qualifier are all handled consistently with the maximum-likelihood path.
#' Left-censored values are substituted at `lod_factor` times the detection limit,
#' right-censored values at the upper limit, and detections kept as reported.
#' Concentrations are converted to mg/kg from the native reporting unit
#' (`ppm` = mg/kg; `pct` = weight percent = 1e4 mg/kg); other units are set to `NA`.
#'
#' @param chem_df Chemistry long table (lab_id, species, unit, qvalue, and
#'   optional qualifier).
#' @param elements Character vector of element symbols to include.
#' @param lod_factor Factor for left-censored substitution (default LOD/2).
#' @param suffix Suffix for column names (default "_mgkg").
#'
#' @return Wide data frame: one row per lab_id, one column per element (mg/kg).
#' @export
build_metals_wide <- function(chem_df,
                              elements,
                              lod_factor = 0.5,
                              suffix     = "_mgkg") {
  if (!"qualifier" %in% names(chem_df)) {
    chem_df <- dplyr::mutate(chem_df, qualifier = NA_character_)
  }

  df <- chem_df %>%
    dplyr::filter(.data$species %in% elements)

  # censoring fields (native units): negative qvalue, .1111 suffix, and L/</G/> qualifier
  cen <- parse_usgs_q(df$qvalue, df$qualifier)

  # substitution value in the native reporting unit
  val_native <- dplyr::case_when(
    cen$left_cen  & is.finite(cen$left_limit)  ~ lod_factor * cen$left_limit,
    cen$right_cen & is.finite(cen$right_limit) ~ cen$right_limit,
    TRUE                                       ~ cen$detected_val
  )

  unit_norm <- stringr::str_to_lower(stringr::str_trim(df$unit))

  # convert native unit -> mg/kg (ppm = mg/kg; pct = weight percent = 1e4 mg/kg)
  to_mgkg <- dplyr::case_when(
    unit_norm == "ppm" ~ 1,
    unit_norm == "pct" ~ 1e4,
    TRUE             ~ NA_real_
  )
  val_mgkg <- val_native * to_mgkg

  tibble::tibble(
    lab_id   = df$lab_id,
    species  = paste0(df$species, suffix),
    val_mgkg = val_mgkg
  ) %>%
    tidyr::pivot_wider(
      id_cols     = "lab_id",
      names_from  = "species",
      values_from = "val_mgkg",
      values_fn   = list(val_mgkg = ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE))
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
      main_df %>% dplyr::select(dplyr::any_of(c("lab_id", "state", "lat_wgs84", "long_wgs84"))),
      by = "lab_id"
    )
  
  if (!is.null(state)) {
    state_code <- norm_state(state)
    out <- out %>% dplyr::filter(.data$state == state_code)
  }
  
  out
}
