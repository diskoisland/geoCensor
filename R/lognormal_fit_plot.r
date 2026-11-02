plot_censored_intervals <- function(df,
                                    col_detected,
                                    col_left_flag,
                                    col_left_limit,
                                    col_right_flag,
                                    col_right_limit,
                                    caller) {
  .censored_intervals(
    df, col_detected, col_left_flag, col_left_limit,
    col_right_flag, col_right_limit, caller
  )
}

#' Plot a censored lognormal fit
#'
#' Fits the same interval-censored lognormal model used by
#' [fit_cen_lognorm()] and plots either the fitted cumulative distribution
#' function or the fitted density on the log10 concentration scale. The
#' geometric mean and 95 percent confidence interval are marked on the
#' concentration axis.
#'
#' @param df Data frame with censoring columns.
#' @param col_detected Column name for detected values.
#' @param col_left_flag Logical column marking left-censored values.
#' @param col_left_limit Numeric LOD for left-censored values.
#' @param col_right_flag Logical column marking right-censored values.
#' @param col_right_limit Numeric upper limit for right-censored values.
#' @param unit_label Label for the concentration unit.
#' @param title Plot title.
#' @param plot_type Plot type: `"cdf"` (default) or `"density"`.
#' @param bins Number of histogram bins for exact detections. Used only when
#'   `plot_type = "density"`.
#' @param show_detects If `TRUE`, show exact detections as a detected-only
#'   empirical CDF scatter in the CDF view.
#' @param show_nonparametric If `TRUE`, show a nonparametric interval-censored
#'   CDF estimate from [survival::survfit()] in the CDF view.
#' @param show_censored If `TRUE`, plot rugs for left- and right-censoring
#'   limits.
#'
#' @return A ggplot object.
#' @export
plot_lognormal_fit <- function(
    df,
    col_detected    = "detected_val",
    col_left_flag   = "left_cen",
    col_left_limit  = "left_limit",
    col_right_flag  = "right_cen",
    col_right_limit = "right_limit",
    unit_label = "ug kg^-1",
    title = "Censored lognormal fit",
    plot_type = c("cdf", "density"),
    bins = 30,
    show_detects = TRUE,
    show_nonparametric = FALSE,
    show_censored = TRUE
) {
  plot_type <- match.arg(plot_type)

  intervals <- plot_censored_intervals(
    df = df,
    col_detected = col_detected,
    col_left_flag = col_left_flag,
    col_left_limit = col_left_limit,
    col_right_flag = col_right_flag,
    col_right_limit = col_right_limit,
    caller = "plot_lognormal_fit"
  )

  if (intervals$n_total < 6 || intervals$n_exact < 2) {
    stop(
      "plot_lognormal_fit: insufficient data; need at least 6 usable ",
      "observations and 2 exact detections"
    )
  }

  fit <- survival::survreg(
    survival::Surv(
      intervals$lower[intervals$ok],
      intervals$upper[intervals$ok],
      type = "interval2"
    ) ~ 1,
    dist = "lognormal"
  )

  meanlog <- unname(stats::coef(fit)[1])
  sdlog <- fit$scale
  se_meanlog <- summary(fit)$table[1, "Std. Error"]
  z <- stats::qnorm(0.975)

  gm <- exp(meanlog)
  ci_low <- exp(meanlog - z * se_meanlog)
  ci_high <- exp(meanlog + z * se_meanlog)

  exact_vals <- intervals$lower[
    intervals$ok &
      is.finite(intervals$lower) &
      intervals$lower > 0 &
      intervals$lower == intervals$upper
  ]
  left_limits <- intervals$upper[
    intervals$ok &
      intervals$lower == 0 &
      is.finite(intervals$upper) &
      intervals$upper > 0
  ]
  right_limits <- intervals$lower[
    intervals$ok &
      is.infinite(intervals$upper) &
      is.finite(intervals$lower) &
      intervals$lower > 0
  ]

  plot_vals <- c(exact_vals, left_limits, right_limits, gm, ci_low, ci_high)
  plot_vals <- plot_vals[is.finite(plot_vals) & plot_vals > 0]
  x_min <- min(plot_vals) / 1.5
  x_max <- max(plot_vals) * 1.5

  meanlog10 <- meanlog / log(10)
  sdlog10 <- sdlog / log(10)

  curve_df <- tibble::tibble(
    log10_value = seq(log10(x_min), log10(x_max), length.out = 400)
  )
  curve_df <- dplyr::mutate(
    curve_df,
    density = stats::dnorm(
      .data$log10_value,
      mean = meanlog10,
      sd = sdlog10
    ),
    cumulative = stats::pnorm(
      .data$log10_value,
      mean = meanlog10,
      sd = sdlog10
    )
  )

  exact_df <- tibble::tibble(value = exact_vals)
  exact_df <- dplyr::mutate(exact_df, log10_value = log10(.data$value))
  exact_cdf_df <- exact_df[order(exact_df$log10_value), , drop = FALSE]
  exact_cdf_df$detected_only_cdf <- stats::ppoints(nrow(exact_cdf_df))

  nonparametric_df <- NULL
  if (identical(plot_type, "cdf") && isTRUE(show_nonparametric)) {
    nonparametric_input <- tibble::tibble(
      lower = intervals$lower[intervals$ok],
      upper = intervals$upper[intervals$ok]
    )
    nonparametric_input <- dplyr::count(
      nonparametric_input,
      .data$lower,
      .data$upper,
      name = "weight"
    )

    nonparametric_fit <- survival::survfit(
      survival::Surv(
        nonparametric_input$lower,
        nonparametric_input$upper,
        type = "interval2"
      ) ~ 1,
      weights = nonparametric_input$weight
    )

    nonparametric_df <- tibble::tibble(
      value = nonparametric_fit$time,
      fitted_cdf = 1 - nonparametric_fit$surv
    )
    nonparametric_df <- dplyr::filter(
      nonparametric_df,
      is.finite(.data$value),
      .data$value > 0,
      is.finite(.data$fitted_cdf)
    )
    nonparametric_df <- dplyr::mutate(
      nonparametric_df,
      log10_value = log10(.data$value)
    )
  }

  fit_label <- sprintf(
    "GM = %.3g %s\n95%% CI: %.3g-%.3g",
    gm,
    unit_label,
    ci_low,
    ci_high
  )

  subtitle <- sprintf(
    "n = %d; exact = %d; left-censored = %d; right-censored = %d",
    intervals$n_total,
    intervals$n_exact,
    intervals$n_left,
    intervals$n_right
  )

  p <- ggplot2::ggplot()

  if (identical(plot_type, "cdf")) {
    if (show_detects && nrow(exact_cdf_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = exact_cdf_df,
          ggplot2::aes(
            x = .data$log10_value,
            y = .data$detected_only_cdf
          ),
          color = "grey35",
          alpha = 0.12,
          size = 0.65,
          inherit.aes = FALSE
        )
    }

    if (!is.null(nonparametric_df) && nrow(nonparametric_df) > 0) {
      p <- p +
        ggplot2::geom_step(
          data = nonparametric_df,
          ggplot2::aes(x = .data$log10_value, y = .data$fitted_cdf),
          color = "grey10",
          linewidth = 0.65,
          alpha = 0.8,
          direction = "hv",
          inherit.aes = FALSE
        )
    }

    p <- p +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$log10_value, y = .data$cumulative),
        color = "#0072B2",
        linewidth = 1
      ) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        color = "grey55",
        linetype = "dotted"
      ) +
      ggplot2::geom_vline(
        xintercept = log10(gm),
        color = "#D55E00",
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = log10(c(ci_low, ci_high)),
        color = "#D55E00",
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::annotate(
        "label",
        x = log10(gm),
        y = 0.08,
        label = fit_label,
        hjust = 0.5,
        size = 3.2,
        fill = "white",
        color = "grey20"
      ) +
      ggplot2::scale_y_continuous(
        name = "Fitted cumulative probability",
        limits = c(0, 1),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = function(x) paste0(round(100 * x), "%")
      ) +
      ggplot2::labs(
        subtitle = subtitle,
        caption = if (show_nonparametric) {
          "Blue: lognormal CDF. Black step: nonparametric censored CDF. Gray: detections only."
        } else {
          "GM = exp(meanlog); gray points are exact detections only."
        }
      )
  } else {
    p <- p +
      ggplot2::geom_histogram(
        data = exact_df,
        ggplot2::aes(
          x = .data$log10_value,
          y = ggplot2::after_stat(.data$density)
        ),
        bins = bins,
        fill = "grey80",
        color = "white",
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = .data$log10_value, y = .data$density),
        color = "#0072B2",
        linewidth = 1
      ) +
      ggplot2::geom_vline(
        xintercept = log10(gm),
        color = "#D55E00",
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = log10(c(ci_low, ci_high)),
        color = "#D55E00",
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::annotate(
        "label",
        x = log10(gm),
        y = Inf,
        label = fit_label,
        vjust = 1.1,
        size = 3.2,
        fill = "white",
        color = "grey20"
      ) +
      ggplot2::labs(
        y = "Density",
        subtitle = subtitle,
        caption = "Histogram shows exact detections only; fitted curve uses exact and censored observations."
      )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      name = paste0("Concentration (", unit_label, ", log10 scale)"),
      labels = function(x) format(signif(10^x, 3), scientific = FALSE)
    ) +
    ggplot2::labs(
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 11)

  if (show_censored) {
    rug_df <- dplyr::bind_rows(
      tibble::tibble(value = left_limits, censoring = "left-censored LOD"),
      tibble::tibble(value = right_limits, censoring = "right-censored limit")
    )
    rug_df <- dplyr::mutate(rug_df, log10_value = log10(.data$value))

    if (nrow(rug_df) > 0) {
      p <- p +
        ggplot2::geom_rug(
          data = rug_df,
          ggplot2::aes(x = .data$log10_value, color = .data$censoring),
          sides = "b",
          alpha = 0.45,
          inherit.aes = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = c(
            "left-censored LOD" = "#009E73",
            "right-censored limit" = "#CC79A7"
          ),
          name = NULL
        )
    }
  }

  p
}
