#' Fit state-level lognormal MLE geomeans (censoring aware)
#'
#' @param ag_main Data frame with state, country and censoring columns
#'   (detected_val, left_cen, left_limit, right_cen, right_limit).
#' @param country_codes Acceptable country values (upper-case).
#'
#' @return Tibble with one row per state: geomean, CI, counts and state-specific
#'   analytical-method review counts/status from [fit_cen_lognorm()].
#' @export
mle_by_state <- function(ag_main,
                         country_codes = c("US", "USA", "UNITED STATES", "U.S.")) {
  ag_main %>%
    dplyr::filter(
      is.na(.data$country) | toupper(.data$country) %in% country_codes
    ) %>%
    dplyr::filter(!is.na(.data$state)) %>%
    dplyr::group_by(.data$state) %>%
    dplyr::group_modify(~ fit_cen_lognorm(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(.data$geomean))
}
