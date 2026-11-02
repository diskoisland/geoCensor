#' Add a substitution column (LOD/2 for left-censored, ULOD for right-censored)
#'
#' @param df Data frame with censoring fields: `detected_val`, `left_cen`,
#'   `left_limit`, `right_cen`, `right_limit` (limits in ug/kg).
#'
#' @return `df` with two added columns: `val_sub_ugkg` and `val_sub_mgkg`,
#'   the substituted value in ug/kg and mg/kg respectively.
#' @export
add_substitution <- function(df) {
  df %>%
    dplyr::mutate(
      val_sub_ugkg = dplyr::case_when(
        .data$left_cen  & is.finite(.data$left_limit)  & .data$left_limit  > 0 ~ 0.5 * .data$left_limit, # LOD/2
        .data$right_cen & is.finite(.data$right_limit) & .data$right_limit > 0 ~ .data$right_limit,       # ULOD
        TRUE ~ .data$detected_val
      ),
      val_sub_mgkg = .data$val_sub_ugkg / 1000
    )
}

#' Censoring-aware lognormal MLE with fitted quantiles
#'
#' Fits an interval-censored lognormal and reports geometric means, model-based
#' confidence limits and fitted concentration quantiles in mg/kg. Input
#' concentrations and limits must be in ug/kg.
#'
#' @param df Data frame with censoring fields: `detected_val`, `left_cen`,
#'   `left_limit`, `right_cen`, `right_limit`.
#' @param probs Unique numeric probabilities between zero and one (inclusive).
#'   Defaults to `c(0.1, 0.5, 0.9, 0.95, 0.99)`. Probability 0.5 is reported as
#'   `median`; other probabilities add `p`-prefixed percentage columns.
#'
#' @details
#' Uses the same interval validation and numerical fit guards as
#' [fit_cen_lognorm()]. Fits require at least six usable records and two exact
#' detections, without requiring any left-censored records. Invalid limits
#' remain excluded; conflicting flags raise an error. Warnings, fitting errors
#' and invalid parameter estimates suppress all fitted values. Count minima
#' are computational guards, not evidence of reliable inference.
#'
#' Confidence limits assume an independent-observation lognormal model; they
#' do not account for spatial dependence or preferential sampling. Spatial
#' confidence limits remain experimental.
#'
#' @return A one-row tibble with `method`, sample counts, `meanlog` (log ug/kg),
#'   `sdlog`, `geomean`, `ci_low`, `ci_high`, `median` and requested quantiles.
#'   Concentration summaries are in mg/kg. The median is always included.
#'   Method is `MLE` on success, `insufficient_data` for insufficient usable
#'   data, or `failed_fit` for a numerical failure. Unavailable fitted values
#'   remain NA, including every requested quantile.
#' @export
fit_lognorm_with_quants <- function(df, probs = c(0.1, 0.5, 0.9, 0.95, 0.99)) {
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || anyDuplicated(probs)) {
    stop("probs must be unique finite numeric probabilities between zero and one")
  }
  quantile_names <- paste0("p", probs[probs != 0.5] * 100)
  if (anyDuplicated(quantile_names)) stop("probs must produce distinct quantile names")
  intervals <- .censored_intervals(df, "detected_val", "left_cen", "left_limit",
    "right_cen", "right_limit", "fit_lognorm_with_quants")
  enough <- intervals$n_total >= 6 && intervals$n_exact >= 2
  fit <- if (enough) .fit_validated_censored_lognormal(
    intervals$lower[intervals$ok], intervals$upper[intervals$ok],
    rep(1, intervals$n_total)
  ) else NULL
  method <- if (!enough) "insufficient_data" else if (is.null(fit)) "failed_fit" else "MLE"
  meanlog <- sdlog <- se_meanlog <- NA_real_
  if (!is.null(fit)) {
    meanlog <- unname(stats::coef(fit)[1L])
    sdlog <- fit$scale
    se_meanlog <- sqrt(stats::vcov(fit)[1L, 1L])
  }
  z <- stats::qnorm(0.975)
  out <- tibble::tibble(
    method = method, n_total = intervals$n_total,
    n_left = intervals$n_left, n_right = intervals$n_right,
    meanlog = meanlog, sdlog = sdlog, geomean = exp(meanlog) / 1000,
    ci_low = exp(meanlog - z * se_meanlog) / 1000,
    ci_high = exp(meanlog + z * se_meanlog) / 1000,
    median = exp(meanlog) / 1000
  )
  values <- if (is.null(fit)) rep(NA_real_, length(quantile_names)) else
    stats::qlnorm(probs[probs != 0.5], meanlog = meanlog, sdlog = sdlog) / 1000
  for (i in seq_along(quantile_names)) out[[quantile_names[i]]] <- values[i]
  # Preserve the established default column order.
  if (identical(probs, c(0.1, 0.5, 0.9, 0.95, 0.99))) {
    out <- out[, c("method", "n_total", "n_left", "n_right", "meanlog", "sdlog",
      "geomean", "ci_low", "ci_high", "p10", "median", "p90", "p95", "p99")]
  }
  out
}
