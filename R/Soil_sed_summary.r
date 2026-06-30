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
#' Fits a left/right interval-censored lognormal (via [survival::survreg()]) and
#' returns the geometric mean, 95% confidence interval, and selected fitted
#' quantiles, all converted from ug/kg to mg/kg.
#'
#' @param df Data frame with censoring fields: `detected_val`, `left_cen`,
#'   `left_limit`, `right_cen`, `right_limit`.
#' @param probs Quantile probabilities to return (default
#'   `c(0.1, 0.5, 0.9, 0.95, 0.99)`).
#'
#' @return A one-row tibble with `geomean`, `ci_low`, `ci_high`, `meanlog`,
#'   `sdlog`, sample counts, and the requested quantiles (mg/kg). Returns NA
#'   values with `method = "MLE"` when there are too few detections to fit.
#' @export
fit_lognorm_with_quants <- function(df, probs = c(0.1, 0.5, 0.9, 0.95, 0.99)) {
  detected_val <- df$detected_val
  left_cen     <- dplyr::coalesce(df$left_cen, FALSE)
  right_cen    <- dplyr::coalesce(df$right_cen, FALSE)
  left_limit   <- df$left_limit
  right_limit  <- df$right_limit

  is_left  <- left_cen  & is.finite(left_limit)  & left_limit  > 0
  is_right <- right_cen & is.finite(right_limit) & right_limit > 0
  is_exact <- is.finite(detected_val) & detected_val > 0 & !is_left & !is_right

  n <- nrow(df)
  lower <- upper <- rep(NA_real_, n)
  lower[is_left]  <- 0;                    upper[is_left]  <- left_limit[is_left]
  lower[is_right] <- right_limit[is_right]; upper[is_right] <- Inf
  lower[is_exact] <- detected_val[is_exact]; upper[is_exact] <- detected_val[is_exact]

  ok <- (is.finite(lower) & is.finite(upper)) | (is.finite(lower) & is.infinite(upper))

  if (sum(ok) < 6 || sum(is_left[ok]) < 2 || sum(is_exact[ok]) < 2) {
    return(tibble::tibble(
      method = "MLE", n_total = sum(ok), n_left = sum(is_left[ok]), n_right = sum(is_right[ok]),
      meanlog = NA_real_, sdlog = NA_real_, geomean = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p10 = NA_real_, median = NA_real_, p90 = NA_real_, p95 = NA_real_, p99 = NA_real_
    ))
  }

  fit     <- survival::survreg(survival::Surv(lower[ok], upper[ok], type = "interval2") ~ 1, dist = "lognormal")
  meanlog <- unname(stats::coef(fit)[1])
  sdlog   <- fit$scale
  qs_ugkg <- stats::qlnorm(probs, meanlog = meanlog, sdlog = sdlog)
  names(qs_ugkg) <- paste0("p", probs * 100)
  z <- stats::qnorm(0.975)

  tibble::tibble(
    method = "MLE", n_total = sum(ok), n_left = sum(is_left[ok]), n_right = sum(is_right[ok]),
    meanlog = meanlog, sdlog = sdlog,
    geomean = exp(meanlog) / 1000,
    ci_low  = exp(meanlog - z * sdlog) / 1000,
    ci_high = exp(meanlog + z * sdlog) / 1000,
    p10 = qs_ugkg["p10"] / 1000, median = qs_ugkg["p50"] / 1000,
    p90 = qs_ugkg["p90"] / 1000, p95 = qs_ugkg["p95"] / 1000, p99 = qs_ugkg["p99"] / 1000
  )
}
