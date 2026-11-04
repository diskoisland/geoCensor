#' Export chemistry records to disk in validated batches
#'
#' Export all records or explicitly selected species without accumulating the
#' selected data in memory. Existing readers and analysis defaults are unchanged.
#'
#' @param chem_csv Path to a readable local chemistry CSV.
#' @param output_dir New output directory whose parent already exists. Existing
#'   files, directories (including empty ones), links and reservations are rejected.
#' @param elements NULL to export every input record, including missing species,
#'   or a unique nonmissing, nonblank character vector of exact species labels.
#'   Matching is case-sensitive and does not trim or harmonize labels.
#' @param chunk_size Maximum parsed records per batch, a positive integer no
#'   larger than `.Machine$integer.max`. Defaults to 10000.
#' @return A list containing `output_dir`, `data_file`, `source_rows_file`,
#'   `info_file`, `files_file` and the import information saved in `info_file`.
#'   No chemistry data frame or vector of all selected positions is returned.
#' @details
#' The new directory contains:
#' \describe{
#'   \item{chemistry.csv}{Selected original columns/records in source order.
#'     Read with [read_usgs_chemistry()] to retain character types. Empty fields
#'     and literal NA follow the reader's missing-marker policy; CSV bytes and
#'     quoting need not match the input. No original column is overwritten.}
#'   \item{source_rows.csv}{An aligned `source_data_record` column of one-based
#'     parsed data-record ordinals, excluding the header and empty lines. Quoted
#'     multiline fields belong to one record. These are not physical line
#'     numbers, sample identifiers or independent-site counts.}
#'   \item{export_info.rds}{Schema/package version, normalized input path,
#'     source byte count and MD5, requested species (NULL means all records),
#'     original column names, input/selected record counts, batch counts and
#'     checksums for the two data files. Use `readRDS()` to inspect it.}
#'   \item{export_files.csv}{Byte counts and MD5 checksums for the two data
#'     files and `export_info.rds`. MD5 is a content fingerprint, not a signature
#'     establishing source authenticity or analytical quality.}
#' }
#' Every batch is parsed as text, including unselected records. Required columns,
#' missing-marker rules and subsequent chemistry-value validation are the same
#' as [read_usgs_chemistry()]. Late malformed records abort without publishing
#' the output directory. Absent species yield header-only data files.
#'
#' A sibling `.NAME.geoCensor-lock` directory reserves the normalized destination
#' for cooperating calls to this function. Data are staged inside this owned
#' reservation and renamed to the new destination only after full parsing,
#' before/after source size/MD5 agreement and provenance writing succeed.
#' Ordinary errors clean up the owned reservation. Forced termination can leave
#' a reservation requiring inspection; existing reservations are never removed
#' automatically. Other programs must leave the source, destination and
#' reservation untouched during export. This is not a crash-durable transaction
#' or protection against noncooperating filesystem mutation.
#'
#' Parser batches are bounded and selected records are written to disk, but
#' total process memory is not guaranteed independent of input size. Output disk
#' space grows with the selected records. No methods, values, QA, dates, sites
#' or dilution factors are inferred or reconciled.
#' @examples
#' path <- tempfile(fileext = ".csv")
#' writeLines(c("lab_id,species,qvalue,unit", "001,Ag,2.0000,ppm",
#'              "001,Pb,10.0000,ppm"), path)
#' destination <- tempfile("chemistry_export_")
#' exported <- export_usgs_chemistry(path, destination, elements = "Ag")
#' read_usgs_chemistry(exported$data_file)
#' readRDS(exported$info_file)$n_input
#' unlink(c(path, destination), recursive = TRUE)
#' @seealso [read_usgs_chemistry()], [read_usgs_chemistry_select()]
#' @export
export_usgs_chemistry <- function(chem_csv, output_dir, elements = NULL, chunk_size = 10000) {
  if (!is.character(chem_csv) || length(chem_csv) != 1L || is.na(chem_csv) ||
      !nzchar(trimws(chem_csv)) || !file.exists(chem_csv) || dir.exists(chem_csv) ||
      file.access(chem_csv, 4L) != 0L) {
    stop("chem_csv must be a readable local CSV file", call. = FALSE)
  }
  if (!is.null(elements) && (!is.character(elements) || !length(elements) ||
      anyNA(elements) || any(!nzchar(trimws(elements))) || anyDuplicated(elements))) {
    stop("elements must be NULL or unique nonmissing, nonblank species labels", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || is.complex(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size) ||
      chunk_size > .Machine$integer.max) {
    stop("chunk_size must be a positive integer no larger than .Machine$integer.max", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(trimws(output_dir))) stop("output_dir must be a new directory path", call. = FALSE)
  output_dir <- path.expand(output_dir)
  leaf <- basename(output_dir)
  if (!nzchar(leaf) || leaf %in% c(".", "..") || !dir.exists(dirname(output_dir))) {
    stop("output_dir must name a new directory with an existing parent", call. = FALSE)
  }
  parent <- normalizePath(dirname(output_dir), winslash = "/", mustWork = TRUE)
  destination <- file.path(parent, leaf)
  if (.export_path_exists(destination)) stop("Output path already exists; choose a new directory", call. = FALSE)
  lock <- file.path(parent, paste0(".", leaf, ".geoCensor-lock"))
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Output destination is reserved or cannot be reserved: ", lock, call. = FALSE)
  }
  token <- basename(tempfile("owner_"))
  owner <- file.path(lock, token)
  on.exit({
    owned <- tryCatch(identical(readLines(owner, warn = FALSE), token), error = function(e) FALSE)
    if (owned) unlink(lock, recursive = TRUE) else {
      warning("Export reservation ownership could not be verified; left for inspection: ", lock, call. = FALSE)
    }
  }, add = TRUE)
  writeLines(token, owner)
  source <- normalizePath(chem_csv, winslash = "/", mustWork = TRUE)
  before <- .export_source_fingerprint(source)
  types <- readr::cols(.default = readr::col_character())
  header <- readr::read_csv(source, n_max = 0, col_types = types,
    na = c("", "NA"), trim_ws = FALSE, name_repair = "check_unique",
    show_col_types = FALSE, progress = FALSE)
  required <- c("lab_id", "species", "qvalue", "unit")
  if (!all(required %in% names(header))) {
    stop("Missing chemistry columns: ", paste(setdiff(required, names(header)), collapse = ", "), call. = FALSE)
  }
  header <- tibble::as_tibble(as.data.frame(header, stringsAsFactors = FALSE))
  stage <- file.path(lock, "contents")
  if (!dir.create(stage)) stop("Could not create export staging directory", call. = FALSE)
  data_file <- file.path(stage, "chemistry.csv")
  rows_file <- file.path(stage, "source_rows.csv")
  readr::write_csv(header, data_file)
  readr::write_csv(tibble::tibble(source_data_record = numeric()), rows_file)
  n_input <- n_selected <- n_chunks <- 0
  max_chunk_rows <- 0L
  callback <- readr::SideEffectChunkCallback$new(function(x, pos) {
    if (!identical(names(x), names(header)) || !all(vapply(x, is.character, logical(1)))) {
      stop("Chemistry CSV chunk schema changed during reading", call. = FALSE)
    }
    if (nrow(readr::problems(x))) {
      stop("Chemistry CSV parsing problems in batch starting at data record ", n_input + 1, call. = FALSE)
    }
    if (!nrow(x)) return(invisible(NULL))
    if (pos != n_input + 1 || nrow(x) > chunk_size) {
      stop("Chemistry CSV parser returned inconsistent record positions or batch size", call. = FALSE)
    }
    keep <- if (is.null(elements)) seq_len(nrow(x)) else which(x$species %in% elements)
    if (length(keep)) {
      readr::write_csv(x[keep, , drop = FALSE], data_file, append = TRUE, col_names = FALSE)
      readr::write_csv(tibble::tibble(source_data_record = n_input + keep), rows_file,
        append = TRUE, col_names = FALSE)
      n_selected <<- n_selected + length(keep)
    }
    n_input <<- n_input + nrow(x)
    n_chunks <<- n_chunks + 1
    max_chunk_rows <<- max(max_chunk_rows, nrow(x))
    invisible(NULL)
  })
  withCallingHandlers(
    readr::read_csv_chunked(source, callback, chunk_size = chunk_size, col_types = types,
      na = c("", "NA"), trim_ws = FALSE, show_col_types = FALSE, progress = FALSE),
    warning = function(w) stop("Chemistry CSV export warning: ", conditionMessage(w), call. = FALSE)
  )
  after <- .export_source_fingerprint(source)
  if (!identical(before, after)) stop("Source file changed during export; no output published", call. = FALSE)
  data_checks <- .export_file_checks(stage, c("chemistry.csv", "source_rows.csv"))
  info <- list(schema_version = 1L, package_version = as.character(getNamespaceVersion("geoCensor")),
    source = list(file = source, bytes = before$bytes, md5 = before$md5),
    elements = elements, columns = names(header), missing_markers = c("", "NA"),
    n_input = n_input, n_selected = n_selected, n_chunks = n_chunks,
    chunk_size = chunk_size, max_chunk_rows = max_chunk_rows, data_files = data_checks)
  saveRDS(info, file.path(stage, "export_info.rds"))
  readr::write_csv(.export_file_checks(stage, c("chemistry.csv", "source_rows.csv", "export_info.rds")),
    file.path(stage, "export_files.csv"))
  if (.export_path_exists(destination)) stop("Output path appeared during export; no output overwritten", call. = FALSE)
  if (!file.rename(stage, destination)) stop("Could not publish export directory", call. = FALSE)
  c(list(output_dir = destination, data_file = file.path(destination, "chemistry.csv"),
    source_rows_file = file.path(destination, "source_rows.csv"),
    info_file = file.path(destination, "export_info.rds"),
    files_file = file.path(destination, "export_files.csv")), info)
}

.export_path_exists <- function(path) {
  link <- suppressWarnings(Sys.readlink(path))
  file.exists(path) || dir.exists(path) || (!is.na(link) && nzchar(link))
}

.export_source_fingerprint <- function(path) {
  bytes <- file.info(path)$size
  md5 <- unname(tools::md5sum(path))
  if (length(bytes) != 1L || is.na(bytes) || is.na(md5)) stop("Could not fingerprint source file", call. = FALSE)
  list(bytes = bytes, md5 = md5)
}

.export_file_checks <- function(directory, files) {
  paths <- file.path(directory, files)
  bytes <- file.info(paths)$size
  md5 <- unname(tools::md5sum(paths))
  if (anyNA(bytes) || anyNA(md5)) stop("Could not fingerprint exported files", call. = FALSE)
  data.frame(file = files, bytes = bytes, md5 = md5, stringsAsFactors = FALSE)
}
