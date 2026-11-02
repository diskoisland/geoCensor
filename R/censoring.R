#' Fit a lognormal with left and right censoring
#'
#' @param df Data frame with censoring columns
#' @param col_detected Column name for detected values
#' @param col_left_flag Logical column marking left-censored
#' @param col_left_limit Numeric LOD for left-censored
#' @param col_right_flag Logical column for right-censored
#' @param col_right_limit Numeric upper limit for right-censored
#'
#' @param weight_col Optional column of nonnegative finite observation weights.
#'   Zero-weight rows are excluded. With weights, confidence limits are omitted:
#'   ordinary model intervals do not account for spatial sampling dependence.
#'
#' @details
#' Fits require at least six usable observations and two exact detections.
#' Uncensored samples and samples with only right censoring are supported.
#' Censored rows with missing, nonpositive, or nonfinite limits are excluded,
#' even when a numeric detected value is present. Rows flagged in both
#' censoring directions raise an error. Missing flags are treated as false.
#' Insufficient samples return `method = "insufficient_data"`; fitting errors,
#' warnings, or invalid parameter estimates return `method = "failed_fit"`.
#' Both cases return missing geometric means and confidence limits.
#'
#' Unweighted confidence limits assume independent observations and an
#' appropriate lognormal model. They do not account for spatial dependence
#' or preferential sampling; spatial confidence limits remain experimental.
#' Minimum observation counts are computational guards, not reliability guarantees.
#'
#' @return A tibble with n_total, n_left, n_right, geomean, ci_low, ci_high, method,
#'   followed by analytical-method review counts and status (see
#'   [flag_analytical_methods()]). These use `amethod` when supplied; absent method
#'   information remains unresolved. Used counts describe eligible fit inputs,
#'   including insufficient or failed fits, and exclude zero-weight rows.
#'
#' @references
#' Helsel, D. R. (1990). Less than obvious: Statistical treatment of data below
#' the detection limit. *Environmental Science & Technology*, 24(12),
#' 1766-1774. \doi{10.1021/es00082a001}
#' @export
fit_cen_lognorm <- function(
    df,
    col_detected    = "detected_val",
    col_left_flag   = "left_cen",
    col_left_limit  = "left_limit",
    col_right_flag  = "right_cen",
    col_right_limit = "right_limit",
    weight_col = NULL
) {
  intervals <- .censored_intervals(
    df, col_detected, col_left_flag, col_left_limit,
    col_right_flag, col_right_limit, "fit_cen_lognorm"
  )
  observation_weights <- rep(1, nrow(df))
  if (!is.null(weight_col)) {
    if (!is.character(weight_col) || length(weight_col) != 1L ||
        is.na(weight_col) || !weight_col %in% names(df)) stop("Missing weight column")
    observation_weights <- df[[weight_col]]
    if (!is.numeric(observation_weights) || any(!is.finite(observation_weights)) ||
        any(observation_weights < 0)) stop("Weights must be numeric, finite, and nonnegative")
    intervals$ok <- intervals$ok & observation_weights > 0
    intervals$n_total <- sum(intervals$ok)
    intervals$n_left <- sum(intervals$ok & intervals$lower == 0, na.rm = TRUE)
    intervals$n_right <- sum(intervals$ok & is.infinite(intervals$upper), na.rm = TRUE)
    intervals$n_exact <- intervals$n_total - intervals$n_left - intervals$n_right
  }
  n_ok <- intervals$n_total
  n_left <- intervals$n_left
  n_right <- intervals$n_right
  n_exact <- intervals$n_exact
  lower <- intervals$lower
  upper <- intervals$upper
  ok <- intervals$ok
  method_summary <- .summarize_analytical_methods(df, intervals = intervals, include_methods = FALSE)$summary
  finish <- function(out) dplyr::bind_cols(out, method_summary)

  if (n_ok < 6 || n_exact < 2) {
    return(finish(tibble::tibble(
      n_total = n_ok,
      n_left  = n_left,
      n_right = n_right,
      geomean = NA_real_,
      ci_low  = NA_real_,
      ci_high = NA_real_,
      method  = "insufficient_data"
    )))
  }
  
  # Normalize to mean one over usable observations; point estimates are invariant
  # to the arbitrary overall scale of the input weights.
  fit_weights <- observation_weights[ok] / max(observation_weights[ok])
  fit_weights <- fit_weights / sum(fit_weights) * n_ok
  fit <- .fit_validated_censored_lognormal(lower[ok], upper[ok], fit_weights)
  if (is.null(fit)) {
    return(finish(tibble::tibble(
      n_total = n_ok, n_left = n_left, n_right = n_right,
      geomean = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      method = "failed_fit"
    )))
  }

  meanlog    <- unname(stats::coef(fit)[1])
  se_meanlog <- sqrt(stats::vcov(fit)[1, 1])
  z          <- stats::qnorm(0.975)
  
  finish(tibble::tibble(
    n_total = n_ok,
    n_left  = n_left,
    n_right = n_right,
    geomean = exp(meanlog),
    ci_low  = if (is.null(weight_col)) exp(meanlog - z * se_meanlog) else NA_real_,
    ci_high = if (is.null(weight_col)) exp(meanlog + z * se_meanlog) else NA_real_,
    method  = if (is.null(weight_col)) "survreg_lognormal_MLE_L+R" else "weighted_survreg_lognormal"
  ))
}


# Shared numerical guard for public lognormal summaries. Censoring intervals
# and minimum observation counts are validated by callers before fitting.
.fit_validated_censored_lognormal <- function(lower, upper, weights) {
  tryCatch(
    withCallingHandlers({
      fit <- survival::survreg(
        survival::Surv(lower, upper, type = "interval2") ~ 1,
        dist = "lognormal", weights = weights
      )
      valid <- all(is.finite(stats::coef(fit))) &&
        is.finite(fit$scale) && fit$scale > 0 &&
        all(is.finite(stats::vcov(fit))) && stats::vcov(fit)[1, 1] > 0
      if (valid) fit else NULL
    }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    error = function(e) NULL
  )
}
