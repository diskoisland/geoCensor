#' Convenience wrapper: read a chemistry CSV and build element censoring fields
#'
#' @param chem_csv Path to a chemistry CSV.
#' @param element Element symbol to filter (default "Ag").
#' @param to_ugkg Logical; if TRUE, convert supported units to ug/kg.
#'
#' @return A data frame retaining rows for the requested element, including
#'   invalid rows with validation_status. Duplicate results are not collapsed.
#' @details All chemistry columns are read as text by [read_usgs_chemistry()],
#'   preserving sparse qualifiers, encoded limits and identifier formatting.
#'   Malformed CSV structure stops before chemistry parsing.
#' @export
build_ag_from_csv <- function(chem_csv, element = "Ag", to_ugkg = TRUE) {
  chem_df <- read_usgs_chemistry(chem_csv)
  build_ag_from_chem(chem_df, element = element, to_ugkg = to_ugkg)
}
