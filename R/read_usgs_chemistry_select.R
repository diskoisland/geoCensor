#' Read selected species from a chemistry CSV in batches
#'
#' Scan the entire file as text, retaining only explicitly requested species and
#' their original data-record positions. The default [read_usgs_chemistry()]
#' continues to read all species at once.
#'
#' @param chem_csv Path to a readable local chemistry CSV.
#' @param elements Unique, nonmissing, nonblank character vector of species
#'   labels. Matching is exact and case-sensitive; input order is retained.
#' @param chunk_size Maximum number of parsed records per batch, a positive
#'   integer no larger than `.Machine$integer.max`. Defaults to 10000.
#' @return A list with:
#' \describe{
#'   \item{data}{A tibble of selected records with all original columns as text,
#'     in source order. Empty fields and the literal NA are missing.}
#'   \item{source_rows}{Numeric vector of one-based parsed data-record ordinals,
#'     aligned with `data`. The header and empty lines are excluded. A quoted
#'     field containing newlines belongs to one record. These are not physical
#'     file line numbers, sample identifiers or independent-site counts.}
#'   \item{n_input}{Total number of parsed input records, including unselected
#'     species.}
#'   \item{n_chunks}{Number of nonempty parser batches scanned.}
#'   \item{max_chunk_rows}{Largest batch scanned, zero for a header-only file.}
#'   \item{elements}{Requested species labels, unchanged.}
#' }
#' @details
#' Required columns and missing-marker rules are the same as
#' [read_usgs_chemistry()]. Every batch is parsed, including unselected records.
#' Duplicate/empty headers, missing required columns and CSV parsing problems
#' raise errors; no partial selection is returned. Chemistry-value validation
#' remains the responsibility of [build_element()] and [parse_usgs_q()].
#'
#' No original column is added, overwritten or renamed; `source_rows` is kept
#' separately even if the source contains a column with that name. Absent species
#' produce a zero-row tibble with the original columns, without observations
#' being invented. Species are not trimmed, harmonized or selected by grade.
#'
#' Parser batches are bounded by `chunk_size`, but retained records and final
#' binding still require memory proportional to the selected output. This is
#' not a constant-memory reader or an assurance that an arbitrary file will fit
#' in memory. The source file must remain unchanged during reading.
#' @examples
#' path <- tempfile(fileext = ".csv")
#' writeLines(c("lab_id,species,qvalue,unit", "001,Ag,2.0000,ppm",
#'              "001,Pb,10.0000,ppm"), path)
#' selected <- read_usgs_chemistry_select(path, "Ag", chunk_size = 1)
#' selected$source_rows
#' build_element(selected$data, element = "Ag")
#' unlink(path)
#' @seealso [read_usgs_chemistry()], [export_usgs_chemistry()], [analyze_elements()]
#' @export
read_usgs_chemistry_select <- function(chem_csv, elements, chunk_size = 10000) {
  if (!is.character(chem_csv) || length(chem_csv) != 1L || is.na(chem_csv) ||
      !nzchar(trimws(chem_csv)) || !file.exists(chem_csv) || dir.exists(chem_csv) ||
      file.access(chem_csv, 4L) != 0L) {
    stop("chem_csv must be a readable local CSV file", call. = FALSE)
  }
  if (!is.character(elements) || !length(elements) || anyNA(elements) ||
      any(!nzchar(trimws(elements))) || anyDuplicated(elements)) {
    stop("elements must contain unique nonmissing, nonblank species labels", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size) ||
      chunk_size > .Machine$integer.max) {
    stop("chunk_size must be a positive integer no larger than .Machine$integer.max", call. = FALSE)
  }
  types <- readr::cols(.default = readr::col_character())
  plain <- function(x) tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
  # The chunked parser has no name_repair argument. Validate before it can
  # repair ambiguous headers, and keep the empty schema for absent species.
  header <- readr::read_csv(chem_csv, n_max = 0, col_types = types,
    na = c("", "NA"), trim_ws = FALSE, name_repair = "check_unique",
    show_col_types = FALSE, progress = FALSE)
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(header))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(header)), collapse = ", "), call. = FALSE)
  }
  header <- plain(header)
  pieces <- positions <- list()
  n_input <- 0
  n_chunks <- 0
  max_chunk_rows <- 0L
  callback <- readr::SideEffectChunkCallback$new(function(x, pos) {
    if (!identical(names(x), names(header)) ||
        !all(vapply(x, is.character, logical(1)))) {
      stop("Chemistry CSV chunk schema changed during reading", call. = FALSE)
    }
    if (nrow(readr::problems(x))) {
      stop("Chemistry CSV parsing problems in batch starting at data record ",
        n_input + 1, ". Resolve malformed records before analysis.", call. = FALSE)
    }
    if (!nrow(x)) return(invisible(NULL))
    if (pos != n_input + 1 || nrow(x) > chunk_size) {
      stop("Chemistry CSV parser returned inconsistent record positions or batch size", call. = FALSE)
    }
    keep <- which(x$species %in% elements)
    if (length(keep)) {
      pieces[[length(pieces) + 1L]] <<- plain(x[keep, , drop = FALSE])
      positions[[length(positions) + 1L]] <<- n_input + keep
    }
    n_input <<- n_input + nrow(x)
    n_chunks <<- n_chunks + 1
    max_chunk_rows <<- max(max_chunk_rows, nrow(x))
    invisible(NULL)
  })
  withCallingHandlers(
    readr::read_csv_chunked(chem_csv, callback, chunk_size = chunk_size,
      col_types = types, na = c("", "NA"), trim_ws = FALSE,
      show_col_types = FALSE, progress = FALSE),
    warning = function(w) stop("Chemistry CSV parsing warning: ", conditionMessage(w), call. = FALSE)
  )
  list(data = if (length(pieces)) dplyr::bind_rows(pieces) else header,
    source_rows = as.numeric(unlist(positions, use.names = FALSE)),
    n_input = n_input, n_chunks = n_chunks, max_chunk_rows = max_chunk_rows,
    elements = elements)
}
