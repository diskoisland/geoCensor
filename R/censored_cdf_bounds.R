#' Bound the empirical CDF without replacing censored values
#'
#' Calculates the minimum and maximum possible empirical cumulative fraction
#' at specified concentrations, retaining left- and right-censoring limits.
#' Optional observation weights can describe an occupied-cell sensitivity
#' analysis. No concentration model is fitted.
#'
#' @param df Data frame with the standard censoring columns.
#' @param thresholds Positive finite concentrations in the same units as df.
#'   Results use the sorted unique thresholds.
#' @param weight_col Optional column of finite nonnegative weights. Zero-weight
#'   rows are excluded. Weights are normalized over usable observations.
#' @param col_detected Column name for exact values.
#' @param col_left_flag,col_left_limit Column names for left-censoring indicators
#'   and upper limits.
#' @param col_right_flag,col_right_limit Column names for right-censoring
#'   indicators and lower limits.
#' @return A list of class geocensor_cdf_bounds, containing cdf (threshold,
#'   cdf_lower, cdf_upper, censoring_gap) and support (usable/excluded counts,
#'   censoring-weight shares and Kish weight effective n). Empty usable subsets
#'   have missing CDF values and weight fractions.
#' @details
#' At a threshold t, exact observations at or below t are counted in both
#' bounds. Left-censored observations count in the lower bound when their limit
#' is at or below t, and can count in the upper bound at every positive t.
#' Right-censored observations never count in the lower bound and can count in
#' the upper bound only when their limit is strictly below t.
#'
#' These are identification bounds on the empirical CDF of the supplied usable
#' records, conditional on the censoring statements and weights. They are not
#' confidence intervals, population CDF bounds, a goodness-of-fit test, or a
#' correction for preferential sampling. Even when both bounds coincide,
#' sampling and spatial uncertainty remain. Weights do not provide concentrations
#' in unobserved cells. Apply within a defined sampling medium/domain.
#' Censored rows follow the same validation/exclusion rules as [fit_cen_lognorm()].
#' @seealso [plot_censored_cdf_bounds()], [declustering_sensitivity()]
#' @export
censored_cdf_bounds <- function(df, thresholds, weight_col = NULL,
    col_detected = "detected_val", col_left_flag = "left_cen",
    col_left_limit = "left_limit", col_right_flag = "right_cen",
    col_right_limit = "right_limit") {
  if (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds <= 0)) {
    stop("thresholds must contain positive finite concentrations")
  }
  thresholds <- sort(unique(thresholds))
  intervals <- .censored_intervals(df, col_detected, col_left_flag,
    col_left_limit, col_right_flag, col_right_limit, "censored_cdf_bounds")
  weights <- rep(1, nrow(df))
  if (!is.null(weight_col)) {
    if (!is.character(weight_col) || length(weight_col) != 1L ||
        is.na(weight_col) || !weight_col %in% names(df)) stop("Missing weight column")
    weights <- df[[weight_col]]
    if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0))
      stop("Weights must be numeric, finite, and nonnegative")
  }
  usable <- intervals$ok & weights > 0
  lower <- intervals$lower[usable]; upper <- intervals$upper[usable]
  exact <- lower == upper; left <- lower == 0; right <- is.infinite(upper)
  n <- sum(usable)
  if (n) {
    weights <- weights[usable] / max(weights[usable])
    weights <- weights / sum(weights)
    # Sorted endpoints avoid an observations-by-thresholds allocation.
    mass_at <- function(values, mass, strict = FALSE) {
      if (!length(values)) return(rep(0, length(thresholds)))
      order <- order(values)
      cumulative <- c(0, cumsum(mass[order]))
      cumulative[findInterval(thresholds, values[order], left.open = strict) + 1L]
    }
    exact_mass <- mass_at(lower[exact], weights[exact])
    minimum <- exact_mass + mass_at(upper[left], weights[left])
    maximum <- exact_mass + sum(weights[left]) +
      mass_at(lower[right], weights[right], strict = TRUE)
    minimum <- pmin(1, pmax(0, minimum))
    maximum <- pmin(1, pmax(minimum, maximum))
    fractions <- c(sum(weights[exact]), sum(weights[left]), sum(weights[right]))
  } else {
    minimum <- maximum <- rep(NA_real_, length(thresholds))
    fractions <- rep(NA_real_, 3L)
    weights <- numeric()
  }
  structure(list(
    cdf = tibble::tibble(threshold = thresholds, cdf_lower = minimum,
      cdf_upper = maximum, censoring_gap = maximum - minimum),
    support = tibble::tibble(n_input = nrow(df), n_total = n,
      n_excluded = sum(!usable), n_exact = sum(exact), n_left = sum(left),
      n_right = sum(right), weighted = !is.null(weight_col),
      exact_fraction = if (n) mean(exact) else NA_real_,
      weighted_exact_fraction = fractions[1L], weighted_left_fraction = fractions[2L],
      weighted_right_fraction = fractions[3L],
      weight_effective_n = if (n) 1 / sum(weights^2) else 0)),
    class = "geocensor_cdf_bounds")
}

#' Plot empirical CDF bounds from censored observations
#'
#' Shows possible empirical cumulative fractions at the evaluated thresholds.
#' Vertical ranges represent unresolved censoring, without imputing values,
#' fitting a distribution or interpolating between thresholds.
#'
#' @param bounds Result from [censored_cdf_bounds()].
#' @param unit_label Concentration unit label.
#' @param title Plot title; use this to identify the sampling medium/domain.
#' @return A ggplot object. No confidence band is calculated.
#' @details
#' This plot describes the supplied records and weights. A narrow censoring
#' range does not establish a precise population estimate. Thresholds are
#' displayed on a log10 concentration axis. Select enough thresholds to inspect
#' the scientifically relevant concentration range; all usable observations
#' contribute to every evaluated bound.
#' @export
plot_censored_cdf_bounds <- function(bounds, unit_label = "ug/kg",
    title = "Empirical CDF censoring bounds") {
  if (!inherits(bounds, "geocensor_cdf_bounds") ||
      !all(c("cdf", "support") %in% names(bounds))) stop("Provide censored_cdf_bounds() output")
  support <- bounds$support
  if (support$n_total < 1L) stop("No usable observations to plot")
  subtitle <- sprintf("Usable = %d; exact = %d; left = %d; right = %d\nExact share of weight = %.2f%%",
    support$n_total, support$n_exact, support$n_left, support$n_right,
    100 * support$weighted_exact_fraction)
  ggplot2::ggplot(bounds$cdf, ggplot2::aes(x = .data$threshold)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$cdf_lower, ymax = .data$cdf_upper),
      color = "grey55", alpha = .45, linewidth = .55) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cdf_lower, color = "Minimum"), size = 1) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cdf_upper, color = "Maximum"), size = 1.3, shape = 1) +
    ggplot2::scale_color_manual(values = c(Minimum = "#2166ac", Maximum = "#b35806"),
      breaks = c("Minimum", "Maximum"), name = "Possible empirical CDF") +
    ggplot2::scale_x_log10(name = paste0("Concentration (", unit_label, ")")) +
    ggplot2::scale_y_continuous(name = "Cumulative fraction of observation weight",
      limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(title = title, subtitle = subtitle,
      caption = "Censoring bounds for the observed records;\nsampling and spatial uncertainty are additional.") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
