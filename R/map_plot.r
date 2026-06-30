#' Ranked barplot of state-level geomeans (log scale)
#'
#' @param mle_state_df Output of mle_by_state().
#' @param unit_label Label to show on y-axis.
#' @param title Plot title.
#'
#' @export
plot_ranked_bar <- function(mle_state_df,
                            unit_label = "µg kg⁻¹",
                            title      = "Silver in U.S. soils by state") {
  ggplot2::ggplot(mle_state_df, ggplot2::aes(
    x = stats::reorder(.data$state, .data$geomean),
    y = .data$geomean
  )) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      width = 0.25
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10(
      name = paste0("Silver (", unit_label, ") — Geomean with 95% CI")
    ) +
    ggplot2::labs(
      x = "State",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' U.S. choropleth map of state-level geomean Ag
#'
#' @param mle_state_df Output of mle_by_state().
#' @param unit_label Label for legend.
#'
#' @export
plot_us_map <- function(mle_state_df,
                        unit_label = "µg kg⁻¹") {
  df <- mle_state_df %>%
    dplyr::mutate(state = toupper(.data$state))
  
  usmap::plot_usmap(data = df, values = "geomean", color = "white") +
    ggplot2::scale_fill_viridis_c(
      name = paste0("Ag (", unit_label, ")"),
      trans = "log10"
    ) +
    ggplot2::labs(
      title = "Geometric mean silver in U.S. soils (state-level MLE)"
    ) +
    ggplot2::theme_minimal()
}

#' State hotspot map for Ag (or any value) in mg/kg
#'
#' @param df Data frame with long_wgs84, lat_wgs84, state, value column.
#' @param state Two-letter code or full name.
#' @param value_col Name of column in df with mg/kg values.
#' @param thr_min Minimum threshold (e.g. 2).
#' @param title Title for the map.
#'
#' @export
state_hotspot_map <- function(df,
                              state,
                              value_col = "Ag_mgkg",
                              thr_min   = 2,
                              title     = NULL) {
  st_code <- norm_state(state)
  if (is.na(st_code)) stop("state_hotspot_map: could not normalize state: ", state)
  
  df_state <- df %>%
    dplyr::filter(.data$state == st_code) %>%
    dplyr::filter(
      is.finite(.data$long_wgs84),
      is.finite(.data$lat_wgs84),
      !!rlang::sym(value_col) >= thr_min
    )
  
  if (nrow(df_state) == 0L) {
    stop("state_hotspot_map: no points ≥ ", thr_min, " mg/kg in ", st_code)
  }
  
  if (is.null(title)) {
    title <- paste0(st_code, " ", value_col, " ≥ ", thr_min, " mg·kg⁻¹")
  }
  
  ggplot2::ggplot(df_state, ggplot2::aes(
    x = .data$long_wgs84,
    y = .data$lat_wgs84,
    color = !!rlang::sym(value_col)
  )) +
    maps::borders("state", regions = tolower(state.name[match(st_code, state.abb)]),
                  colour = "grey60", fill = NA) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::coord_quickmap() +
    ggplot2::scale_color_viridis_c(
      name  = paste0("Ag (mg kg", "\u2212", "1)"),
      trans = "log10"
    ) +
    ggplot2::labs(
      title = title,
      x = "Longitude", y = "Latitude"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}