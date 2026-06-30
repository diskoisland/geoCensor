#' Convenience wrapper: read a chemistry CSV and build Ag censoring fields
#' @export
build_ag_from_csv <- function(chem_csv, element = "Ag", to_ugkg = TRUE) {
  chem_df <- readr::read_csv(chem_csv, show_col_types = FALSE)
  build_ag_from_chem(chem_df, element = element, to_ugkg = to_ugkg)
}