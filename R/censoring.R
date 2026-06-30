#' Fit a lognormal with left and right censoring
#'
#' @param df Data frame with censoring columns
#' @param col_detected Column name for detected values
#' @param col_left_flag Logical column marking left-censored
#' @param col_left_limit Numeric LOD for left-censored
#' @param col_right_flag Logical column for right-censored
#' @param col_right_limit Numeric upper limit for right-censored
#'
#' @return A tibble with n_total, n_left, n_right, geomean, ci_low, ci_high, method.
#' @export
fit_cen_lognorm <- function(
    df,
    col_detected    = "detected_val",
    col_left_flag   = "left_cen",
    col_left_limit  = "left_limit",
    col_right_flag  = "right_cen",
    col_right_limit = "right_limit"
) {
  req <- c(col_detected, col_left_flag, col_left_limit,
           col_right_flag, col_right_limit)
  if (!all(req %in% names(df))) {
    stop("fit_cen_lognorm: Missing columns: ",
         paste(req[!req %in% names(df)], collapse = ", "))
  }
  
  detected_val <- df[[col_detected]]
  left_flag    <- df[[col_left_flag]]
  right_flag   <- df[[col_right_flag]]
  left_limit   <- df[[col_left_limit]]
  right_limit  <- df[[col_right_limit]]
  
  left_cen  <- as.logical(left_flag);  left_cen[is.na(left_cen)]   <- FALSE
  right_cen <- as.logical(right_flag); right_cen[is.na(right_cen)] <- FALSE
  
  is_left  <- left_cen  & is.finite(left_limit)  & left_limit  > 0
  is_right <- right_cen & is.finite(right_limit) & right_limit > 0
  is_exact <- is.finite(detected_val) & detected_val > 0 &
    !is_left & !is_right
  
  n <- nrow(df)
  lower <- upper <- rep(NA_real_, n)
  
  lower[is_left]  <- 0
  upper[is_left]  <- left_limit[is_left]
  
  lower[is_right] <- right_limit[is_right]
  upper[is_right] <- Inf
  
  lower[is_exact] <- detected_val[is_exact]
  upper[is_exact] <- detected_val[is_exact]
  
  ok      <- (is_left | is_right | is_exact)
  n_ok    <- sum(ok)
  n_left  <- sum(is_left)
  n_right <- sum(is_right)
  n_exact <- sum(is_exact)
  
  if (n_ok < 6 || n_left < 1 || n_exact < 2) {
    return(tibble::tibble(
      n_total = n_ok,
      n_left  = n_left,
      n_right = n_right,
      geomean = NA_real_,
      ci_low  = NA_real_,
      ci_high = NA_real_,
      method  = "insufficient_data"
    ))
  }
  
  fit <- survival::survreg(
    survival::Surv(lower[ok], upper[ok], type = "interval2") ~ 1,
    dist = "lognormal"
  )
  
  meanlog    <- unname(stats::coef(fit)[1])
  se_meanlog <- summary(fit)$table[1, "Std. Error"]
  z          <- stats::qnorm(0.975)
  
  tibble::tibble(
    n_total = n_ok,
    n_left  = n_left,
    n_right = n_right,
    geomean = exp(meanlog),
    ci_low  = exp(meanlog - z * se_meanlog),
    ci_high = exp(meanlog + z * se_meanlog),
    method  = "survreg_lognormal_MLE_L+R"
  )
}