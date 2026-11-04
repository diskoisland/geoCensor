#' Build a wide table of metals (mg/kg) using substitution for censoring
#'
#' Censoring is parsed with [parse_usgs_q()], so left/right censoring encoded via
#' negative qvalues, the `.1111` upper-limit suffix, **or** the `L`/`<`/`G`/`>`
#' qualifier are all handled consistently with the maximum-likelihood path.
#' Left-censored values are substituted at `lod_factor` times the detection limit,
#' right-censored values at the upper limit, and detections kept as reported.
#' Concentrations are converted to mg/kg from the native reporting unit
#' (ppm, mg/kg, ppb, ug/kg, micro-symbol equivalents, pct, percent, or %).
#' Other units are set to `NA`. Duplicate sample/element pairs raise an error.
#' Rounded summary upper-limit recovery uses the same known-method and coding
#' guard as [build_element()]. Inspect its long observations for row-level
#' `summary_recovery_status` before widening.
#'
#' @param chem_df Chemistry long table (lab_id, species, unit, qvalue, and
#'   optional qualifier).
#' @param elements Character vector of element symbols to include.
#' @param lod_factor Factor for left-censored substitution (default LOD/2).
#' @param suffix Suffix for column names (default "_mgkg").
#'
#' @return Wide data frame: one row per lab_id, one column per element (mg/kg),
#'   with a `validation_summary` attribute counting usable and excluded rows.
#'   Excluded input rows produce a warning and missing concentrations.
#' @export
build_metals_wide <- function(chem_df,
                              elements,
                              lod_factor = 0.5,
                              suffix     = "_mgkg") {
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(chem_df))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(chem_df)), collapse = ", "))
  }
  if (!is.numeric(lod_factor) || length(lod_factor) != 1L ||
      !is.finite(lod_factor) || lod_factor <= 0 || lod_factor > 1) {
    stop("lod_factor must be a number greater than zero and at most one")
  }
  if (!"qualifier" %in% names(chem_df)) {
    chem_df <- dplyr::mutate(chem_df, qualifier = NA_character_)
  }

  df <- chem_df %>%
    dplyr::filter(.data$species %in% elements)

  # censoring fields (native units): negative qvalue, .1111 suffix, and L/</G/> qualifier
  cen <- parse_usgs_q(.qvalue_for_censoring(df)$value, df$qualifier)

  # substitution value in the native reporting unit
  val_native <- dplyr::case_when(
    cen$left_cen  & is.finite(cen$left_limit)  ~ lod_factor * cen$left_limit,
    cen$right_cen & is.finite(cen$right_limit) ~ cen$right_limit,
    TRUE                                       ~ cen$detected_val
  )

  if (anyDuplicated(df[c("lab_id", "species")])) {
    stop("Duplicate lab_id/species rows; select one result per sample and element before widening")
  }
  to_mgkg <- .concentration_factor_ugkg(df$unit) / 1000
  val_mgkg <- val_native * to_mgkg

  validation_status <- cen$validation_status
  validation_status[validation_status == "usable" & !is.finite(to_mgkg)] <- "unsupported_unit"
  if (any(validation_status != "usable")) {
    warning(sum(validation_status != "usable"),
            " chemistry rows excluded; see the validation_summary attribute", call. = FALSE)
  }
  out <- tibble::tibble(
    lab_id   = df$lab_id,
    species  = paste0(df$species, suffix),
    val_mgkg = val_mgkg
  ) %>%
    tidyr::pivot_wider(
      id_cols     = "lab_id",
      names_from  = "species",
      values_from = "val_mgkg"
    )
  attr(out, "validation_summary") <- dplyr::count(
    tibble::tibble(validation_status = validation_status),
    .data$validation_status, name = "n_observations"
  )
  out
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
  
  out <- join_soil_main(
    wide, main_df, extra_cols = c("state", "lat_wgs84", "long_wgs84")
  )

  if (!is.null(state)) {
    state_code <- norm_state(state)
    out <- out %>% dplyr::filter(.data$state == state_code)
  }
  
  out
}
