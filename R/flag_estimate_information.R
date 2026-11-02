#' Flag count-based information limits in concentration estimates
#'
#' Adds exact-observation counts and configurable review flags without changing
#' any estimate or excluding rows. These are descriptive review thresholds,
#' not validated requirements for reliable inference.
#'
#' @param df Estimate table with n_total, n_left, n_right, and geomean.
#' @param min_exact Minimum exact observations before the few-exact review flag
#'   clears. Default 10 is a review convention, not a statistical guarantee.
#' @param min_exact_fraction Minimum exact-observation fraction before the
#'   low-fraction review flag clears. Default 0.05 is configurable.
#' @return Input table with n_exact, exact_fraction, review_few_exact,
#'   review_low_exact_fraction, and review_status. Status is estimate_unavailable,
#'   review_information, or counts_above_review_thresholds. Existing derived
#'   columns with these names are recalculated. Estimates are unchanged.
#' @details
#' Exact observations exclude both left and right censoring; right-censored
#' observations are not being classified as nondetections. A zero usable count
#' gives a missing exact fraction. Count flags do not assess fit adequacy,
#' spatial independence, support comparability, preferential selection, or bias.
#' Clearing the flags must not be interpreted as a certification of validity.
#' @export
flag_estimate_information <- function(df, min_exact = 10, min_exact_fraction = 0.05) {
  required <- c("n_total", "n_left", "n_right", "geomean")
  if (!all(required %in% names(df))) stop("Missing estimate count or geometric mean columns")
  if (!is.numeric(min_exact) || length(min_exact) != 1L || !is.finite(min_exact) ||
      min_exact < 0 || min_exact != floor(min_exact)) stop("min_exact must be one nonnegative whole number")
  if (!is.numeric(min_exact_fraction) || length(min_exact_fraction) != 1L ||
      !is.finite(min_exact_fraction) || min_exact_fraction < 0 || min_exact_fraction > 1) {
    stop("min_exact_fraction must be in [0, 1]")
  }
  for (column in c("n_total", "n_left", "n_right")) {
    value <- df[[column]]
    if (!is.numeric(value) || any(!is.finite(value)) || any(value < 0) || any(value != floor(value))) {
      stop("Observation counts must be finite nonnegative whole numbers")
    }
  }
  if (any(df$n_left + df$n_right > df$n_total)) stop("Censored counts exceed n_total")
  if (!is.numeric(df$geomean)) stop("geomean must be numeric")
  out <- tibble::as_tibble(df)
  out$n_exact <- out$n_total - out$n_left - out$n_right
  out$exact_fraction <- ifelse(out$n_total > 0, out$n_exact / out$n_total, NA_real_)
  out$review_few_exact <- out$n_exact < min_exact
  out$review_low_exact_fraction <- out$exact_fraction < min_exact_fraction
  unavailable <- !is.finite(out$geomean) | out$geomean <= 0
  review <- out$review_few_exact | (!is.na(out$review_low_exact_fraction) & out$review_low_exact_fraction)
  out$review_status <- ifelse(unavailable, "estimate_unavailable",
    ifelse(review, "review_information", "counts_above_review_thresholds"))
  out
}
