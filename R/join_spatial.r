#' Join element data (e.g., Ag) to NGDB main metadata by lab_id
#'
#' @param elem_df Data frame with at least lab_id.
#' @param main_df Output of \code{read_clean_main()}.
#' @param extra_cols Character vector of extra metadata cols to keep.
#'
#' @return Joined data frame with a `join_audit` attribute containing matched,
#'   unmatched, and repeated chemistry row counts. Duplicate or missing metadata
#'   IDs raise an error. Unmatched and repeated chemistry rows warn and are retained.
#' @export
join_soil_main <- function(elem_df,
                           main_df,
                           extra_cols = c(
                             "state", "country",
                             "lat_wgs84", "long_wgs84",
                             "depth_num", "depth_unit"
                           )) {
  for (x in list(elem_df, main_df)) {
    if (!"lab_id" %in% names(x)) stop("Both tables must contain lab_id")
  }
  elem_df$lab_id <- trimws(as.character(elem_df$lab_id))
  main_df$lab_id <- trimws(as.character(main_df$lab_id))
  missing_main <- is.na(main_df$lab_id) | main_df$lab_id == ""
  if (any(missing_main)) stop("Metadata lab_id values must not be missing or blank")
  if (anyDuplicated(main_df$lab_id)) {
    stop("Duplicate metadata lab_id values would multiply observations; resolve duplicates before joining")
  }
  unmatched <- is.na(elem_df$lab_id) | elem_df$lab_id == "" |
    !elem_df$lab_id %in% main_df$lab_id
  if (any(unmatched)) {
    warning(sum(unmatched), " chemistry rows have no matching metadata", call. = FALSE)
  }
  duplicate_rows <- sum(duplicated(elem_df$lab_id[!is.na(elem_df$lab_id) & elem_df$lab_id != ""]))
  if (duplicate_rows > 0L) {
    warning(duplicate_rows, " repeated chemistry lab_id rows retained; verify sample independence", call. = FALSE)
  }
  main_keep <- main_df %>%
    dplyr::select(dplyr::any_of(c("lab_id", extra_cols)))
  out <- dplyr::left_join(elem_df, main_keep, by = "lab_id", na_matches = "never")
  attr(out, "join_audit") <- tibble::tibble(
    n_chemistry = nrow(elem_df), n_matched = sum(!unmatched),
    n_unmatched = sum(unmatched), n_repeated_chemistry = duplicate_rows
  )
  out
}

#' Alias for joining sediments (identical behavior)
#'
#' @param elem_df Data frame with at least lab_id.
#' @param main_df Output of \code{read_clean_main()}.
#' @param extra_cols Character vector of extra metadata cols to keep.
#'
#' @return Joined data frame.
#' @export
join_sed_main <- function(elem_df,
                          main_df,
                          extra_cols = c(
                            "state", "country",
                            "lat_wgs84", "long_wgs84",
                            "depth_num", "depth_unit"
                          )) {
  join_soil_main(elem_df, main_df, extra_cols = extra_cols)
}
