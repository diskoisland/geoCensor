#' Read and clean USGS NGDB main metadata
#'
#' @param main_csv Path to main.csv (soil or sediment).
#'
#' @return A cleaned data frame with normalized state, numeric lat/long, and depth.
#' @export
read_clean_main <- function(main_csv) {
  ct <- readr::cols(
    lab_id              = readr::col_character(),
    field_id            = readr::col_character(),
    job_id              = readr::col_character(),
    submitter           = readr::col_character(),
    date_collect        = readr::col_character(),
    date_submitted      = readr::col_character(),
    project_name        = readr::col_character(),
    state               = readr::col_character(),
    country             = readr::col_character(),
    lat_wgs84           = readr::col_double(),
    long_wgs84          = readr::col_double(),
    datum_orig          = readr::col_character(),
    lat_orig            = readr::col_double(),
    long_orig           = readr::col_double(),
    coord_prec_est      = readr::col_character(),
    depth               = readr::col_character(),
    locate_desc         = readr::col_character(),
    sample_source       = readr::col_character(),
    method_collected    = readr::col_character(),
    primary_class       = readr::col_character(),
    specific_name       = readr::col_character(),
    addl_attr           = readr::col_character(),
    geologic_age        = readr::col_character(),
    stratigraphy        = readr::col_character(),
    mineralization      = readr::col_character(),
    prep                = readr::col_character(),
    mesh_size           = readr::col_character(),
    collection_protocol = readr::col_character(),
    previous_job_id     = readr::col_character(),
    previous_lab_id     = readr::col_character()
  )
  
  main <- readr::read_csv(
    main_csv,
    col_types       = ct,
    na              = c("", "NA", "NaN"),
    show_col_types  = FALSE
  )
  
  main %>%
    dplyr::mutate(
      state   = norm_state(.data$state),
      country = stringr::str_trim(.data$country),
      
      lat_wgs84 = dplyr::if_else(
        is.finite(.data$lat_wgs84) & .data$lat_wgs84 >= -90 & .data$lat_wgs84 <= 90,
        .data$lat_wgs84, NA_real_
      ),
      long_wgs84 = dplyr::if_else(
        is.finite(.data$long_wgs84) & .data$long_wgs84 >= -180 & .data$long_wgs84 <= 180,
        .data$long_wgs84, NA_real_
      ),
      
      date_collect   = suppressWarnings(lubridate::mdy(.data$date_collect)),
      date_submitted = suppressWarnings(lubridate::mdy(.data$date_submitted)),
      
      depth_num = {
        rng_hi <- suppressWarnings(
          readr::parse_number(stringr::str_extract(.data$depth, "(?<=-)[0-9.]+"))
        )
        one <- suppressWarnings(readr::parse_number(.data$depth))
        dplyr::coalesce(rng_hi, one)
      },
      depth_unit = dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(.data$depth), "cm") ~ "cm",
        stringr::str_detect(stringr::str_to_lower(.data$depth), "mm") ~ "mm",
        stringr::str_detect(stringr::str_to_lower(.data$depth), "\\bm\\b") ~ "m",
        TRUE ~ NA_character_
      )
    )
}