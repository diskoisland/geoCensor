#' Read chemistry CSV fields without guessing their types
#'
#' Preserve identifiers, reported values, qualifiers, methods and source-summary
#' text before passing them to the censoring-aware builders.
#'
#' @param chem_csv Path to a readable local chemistry CSV.
#' @return A tibble with every column read as character, in original row/column
#'   order. Empty fields and the literal NA are missing; other field text,
#'   including whitespace and lexical numeric precision, is retained.
#' @details
#' Required columns are lab_id, species, qvalue and unit. Duplicate/empty headers,
#' malformed CSV records and missing required columns raise errors. Individual
#' chemistry values and qualifiers are validated later by [build_element()] and
#' [parse_usgs_q()]; invalid values are not guessed, repaired or replaced here.
#'
#' Reading all columns as text prevents sparse qualifiers or encoded upper limits
#' from being lost through numeric/logical type guessing. A positive qvalue with
#' a rare G qualifier must remain an upper limit, even late in a large file.
#' Leading-zero and large identifiers are text, not numbers. The reader does not
#' select results, reconcile methods, infer dates or certify analytical quality.
#' @seealso [read_usgs_chemistry_select()], [build_ag_from_csv()], [analyze_elements()]
#' @export
read_usgs_chemistry <- function(chem_csv) {
  if (!is.character(chem_csv) || length(chem_csv) != 1L || is.na(chem_csv) ||
      !nzchar(trimws(chem_csv)) || !file.exists(chem_csv) || dir.exists(chem_csv)) {
    stop("chem_csv must be a readable local CSV file")
  }
  data <- readr::read_csv(chem_csv,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA"), trim_ws = FALSE, name_repair = "check_unique",
    show_col_types = FALSE, progress = FALSE)
  problems <- readr::problems(data)
  if (nrow(problems)) {
    stop("Chemistry CSV parsing problems at file rows: ",
      paste(utils::head(unique(problems$row), 10L), collapse = ", "),
      ". Resolve malformed records before analysis.", call. = FALSE)
  }
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(data))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  tibble::as_tibble(as.data.frame(data, stringsAsFactors = FALSE))
}
