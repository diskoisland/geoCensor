#' Join element data (e.g., Ag) to NGDB main metadata by lab_id
#'
#' @param elem_df Data frame with at least lab_id.
#' @param main_df Output of \code{read_clean_main()}.
#' @param extra_cols Character vector of extra metadata cols to keep.
#'
#' @return Joined data frame.
#' @export
join_soil_main <- function(elem_df,
                           main_df,
                           extra_cols = c(
                             "state", "country",
                             "lat_wgs84", "long_wgs84",
                             "depth_num", "depth_unit"
                           )) {
  main_keep <- main_df %>%
    dplyr::select(dplyr::any_of(c("lab_id", extra_cols)))
  
  elem_df %>%
    dplyr::left_join(main_keep, by = "lab_id")
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
