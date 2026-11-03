#' Save a reproducible sampling-support audit from a metadata CSV
#'
#' Read metadata as text, audit recorded support and optionally export exact-position
#' memberships. Preserve source field text separately from numeric coordinates.
#'
#' @param input_file Local metadata CSV path.
#' @param output_dir New or empty destination directory.
#' @param group_cols,support_cols,id_col Passed to [sampling_support_audit()].
#' @param x_col,y_col Optional coordinate columns. Both are required together.
#' @param include_memberships Append record memberships and global coordinate
#'   groups; requires coordinates. Default FALSE.
#' @return Invisibly, a list with input_metadata (raw text), audit (the support
#'   audit tables), and coordinate_conversion (raw/numeric coordinate provenance).
#'   The same object is saved as audit.rds alongside CSV tables, run_summary.txt
#'   and session_info.txt.
#' @details
#' CSV fields are read as character without whitespace trimming or type guessing.
#' Raw field values and leading-zero identifiers are retained in input_metadata
#' and optional memberships. In the analysis copy, whitespace-only fields and
#' the trimmed literal NA are missing. Other labels are not reconciled.
#'
#' Coordinates accept decimal/scientific numeric text, blank/NA missing values,
#' and explicit Inf/-Inf/NaN. Invalid text aborts before output creation; no CRS,
#' unit or degree-minute conversion is inferred. Nonfinite values remain unlocated.
#' Coordinate conversion is documented per input row/column. Membership source
#' columns retain original text; coordinate-group/location columns are numeric.
#'
#' Output includes input_metadata.csv, fields.csv, combinations.csv and, when
#' coordinates are supplied, locations.csv and coordinate_conversion.csv. Optional
#' memberships.csv and coordinate_groups.csv appear only when requested. Input
#' MD5 checks before/after analysis detect source changes before saving. Existing
#' nonempty destinations are rejected. These are descriptive review outputs,
#' not physical-site identification, comparability or independence certification.
#' @export
audit_sampling_file <- function(input_file, output_dir, group_cols = character(),
                                support_cols = NULL, id_col = "lab_id",
                                x_col = NULL, y_col = NULL,
                                include_memberships = FALSE) {
  path_ok <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  if (!path_ok(input_file) || !file.exists(input_file) || dir.exists(input_file)) stop("input_file must be a readable local CSV file")
  if (!path_ok(output_dir)) stop("output_dir must name one directory")
  if (!is.logical(include_memberships) || length(include_memberships) != 1L || is.na(include_memberships)) {
    stop("include_memberships must be TRUE or FALSE")
  }
  if (xor(is.null(x_col), is.null(y_col))) stop("Supply both coordinate columns")
  if (include_memberships && is.null(x_col)) stop("include_memberships requires both coordinate columns")
  if (!is.null(x_col) && (!path_ok(x_col) || !path_ok(y_col) || x_col == y_col)) stop("Coordinates must name two distinct columns")
  check_destination <- function() {
    if (file.exists(output_dir) && !dir.exists(output_dir)) stop("Output path is not a directory")
    if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
      stop("Output directory is not empty; choose a new directory to preserve previous results")
    }
  }
  check_destination()
  checksum <- unname(tools::md5sum(input_file))
  if (is.na(checksum)) stop("Could not fingerprint input_file")
  raw <- readr::read_csv(input_file, col_types = readr::cols(.default = readr::col_character()),
    na = character(), trim_ws = FALSE, name_repair = "check_unique", show_col_types = FALSE, progress = FALSE)
  if (nrow(readr::problems(raw))) stop("CSV parsing problems; resolve malformed input before auditing")
  raw <- tibble::as_tibble(as.data.frame(raw, stringsAsFactors = FALSE))
  if (!path_ok(id_col) || !id_col %in% names(raw)) stop("Missing sample identifier column")
  data <- raw
  for (column in names(data)) {
    missing <- is.na(data[[column]]) | trimws(data[[column]]) %in% c("", "NA")
    data[[column]][missing] <- NA_character_
  }
  conversions <- tibble::tibble(input_row = integer(), identifier_raw = character(),
    coordinate_column = character(), raw_value = character(), numeric_value = double(),
    conversion_status = character())
  if (!is.null(x_col)) {
    if (!all(c(x_col, y_col) %in% names(raw))) stop("Missing coordinate columns")
    decimal <- "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
    for (column in c(x_col, y_col)) {
      text <- trimws(raw[[column]])
      blank <- is.na(text) | text == ""
      missing_na <- !is.na(text) & text == "NA"
      special <- !is.na(text) & toupper(text) %in% c("INF", "+INF", "-INF", "NAN")
      numeric_text <- !is.na(text) & grepl(decimal, text, perl = TRUE)
      invalid <- which(!(blank | missing_na | special | numeric_text))
      if (length(invalid)) stop("Coordinate column '", column, "' has nonnumeric values at input rows: ",
        paste(utils::head(invalid, 10L), collapse = ", "))
      value <- rep(NA_real_, length(text))
      value[numeric_text] <- suppressWarnings(as.numeric(text[numeric_text]))
      lookup <- c(INF = Inf, "+INF" = Inf, "-INF" = -Inf, NAN = NaN)
      value[special] <- unname(lookup[toupper(text[special])])
      status <- rep("nonfinite", length(text))
      status[is.finite(value)] <- "parsed_finite"
      status[blank] <- "missing_blank"
      status[missing_na] <- "missing_NA"
      data[[column]] <- value
      conversions <- dplyr::bind_rows(conversions, tibble::tibble(
        input_row = seq_len(nrow(raw)), identifier_raw = raw[[id_col]], coordinate_column = column,
        raw_value = raw[[column]], numeric_value = value, conversion_status = status))
    }
  }
  audit <- sampling_support_audit(data, group_cols = group_cols, support_cols = support_cols,
    id_col = id_col, x_col = x_col, y_col = y_col, include_memberships = include_memberships)
  if (include_memberships) {
    for (column in names(raw)) audit$memberships[[column]] <- raw[[column]]
  }
  if (!identical(unname(tools::md5sum(input_file)), checksum)) stop("Input file changed during analysis; no output saved")
  check_destination()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create output directory")
  readr::write_csv(raw, file.path(output_dir, "input_metadata.csv"), na = "")
  for (name in c("fields", "combinations", if (!is.null(x_col)) "locations",
                 if (include_memberships) c("memberships", "coordinate_groups"))) {
    readr::write_csv(audit[[name]], file.path(output_dir, paste0(name, ".csv")), na = "")
  }
  if (!is.null(x_col)) readr::write_csv(conversions, file.path(output_dir, "coordinate_conversion.csv"), na = "")
  result <- list(input_metadata = raw, audit = audit, coordinate_conversion = conversions)
  saveRDS(result, file.path(output_dir, "audit.rds"))
  writeLines(c(paste("Input:", normalizePath(input_file)), paste("Input MD5:", checksum),
    paste("Package version:", utils::packageVersion("geoCensor")),
    paste("Records:", nrow(raw)), paste("Identifier column:", id_col),
    paste("Analysis groups:", if (length(group_cols)) paste(group_cols, collapse = ", ") else "pooled"),
    paste("Support fields:", paste(audit$available_fields, collapse = ", ")),
    paste("Coordinates:", if (is.null(x_col)) "not supplied" else paste(x_col, y_col, sep = ", ")),
    paste("Memberships:", include_memberships),
    if (!is.null(x_col)) paste("Unlocated records:", sum(!is.finite(data[[x_col]]) | !is.finite(data[[y_col]]))),
    "Analysis missing markers: trimmed blank and literal NA; raw input text is retained.",
    "Coordinate system and units are supplied by the user; no transformation is inferred.",
    "Shared positions are review groups, not certified physical sites, replicates or independent observations.",
    "No records collapsed, weights assigned or dilution corrections applied."), file.path(output_dir, "run_summary.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "session_info.txt"))
  invisible(result)
}
