#' Normalize state names to USPS codes
#'
#' @param x Character vector of state names or abbreviations
#'
#' @return Uppercase 2-letter state abbreviations, or NA
#' @export
norm_state <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  
  state_map <- stats::setNames(datasets::state.abb, datasets::state.name)
  
  is_two <- !is.na(x) & nchar(x) == 2
  x[is_two] <- toupper(x[is_two])
  
  is_name <- !is_two & !is.na(x)
  if (any(is_name)) {
    full <- stringr::str_to_title(x[is_name])
    mapped <- state_map[full]
    x[is_name] <- mapped
  }
  
  x[!x %in% datasets::state.abb] <- NA_character_
  x
}
