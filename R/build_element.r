#' Build censoring-aware data for a single element (generic)
#'
#' @param chem_df Chemistry data frame.
#' @param element Element symbol (e.g. "Pb", "Zn").
#' @param to_ugkg Convert ppm → µg/kg?
#'
#' @return Data frame; same structure as \code{build_ag_from_chem()} but for any element.
#' @export
build_element <- function(chem_df,
                          element,
                          to_ugkg = TRUE) {
  build_ag_from_chem(chem_df, element = element, to_ugkg = to_ugkg)
}