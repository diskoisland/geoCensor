.censored_intervals <- function(df,
                                col_detected,
                                col_left_flag,
                                col_left_limit,
                                col_right_flag,
                                col_right_limit,
                                caller) {
  req <- c(col_detected, col_left_flag, col_left_limit,
           col_right_flag, col_right_limit)
  if (!all(req %in% names(df))) {
    stop(caller, ": Missing columns: ",
         paste(req[!req %in% names(df)], collapse = ", "))
  }

  detected_val <- df[[col_detected]]
  left_flag    <- df[[col_left_flag]]
  right_flag   <- df[[col_right_flag]]
  left_limit   <- df[[col_left_limit]]
  right_limit  <- df[[col_right_limit]]

  left_cen  <- as.logical(left_flag)
  right_cen <- as.logical(right_flag)
  left_cen[is.na(left_cen)] <- FALSE
  right_cen[is.na(right_cen)] <- FALSE

  is_left  <- left_cen & is.finite(left_limit) & left_limit > 0
  is_right <- right_cen & is.finite(right_limit) & right_limit > 0
  is_exact <- is.finite(detected_val) & detected_val > 0 & !is_left & !is_right

  n <- nrow(df)
  lower <- upper <- rep(NA_real_, n)

  lower[is_left]  <- 0
  upper[is_left]  <- left_limit[is_left]

  lower[is_right] <- right_limit[is_right]
  upper[is_right] <- Inf

  lower[is_exact] <- detected_val[is_exact]
  upper[is_exact] <- detected_val[is_exact]

  ok <- is_left | is_right | is_exact

  list(
    lower = lower,
    upper = upper,
    ok = ok,
    n_total = sum(ok),
    n_left = sum(is_left[ok]),
    n_right = sum(is_right[ok]),
    n_exact = sum(is_exact[ok])
  )
}

#' Compare censored parametric distributions
#'
#' Fits one or more parametric distributions to left/right/interval-censored
#' concentration data using [survival::survreg()] and compares them by AIC.
#' This is a diagnostic for whether the lognormal assumption is competitive
#' against other censored parametric choices such as Weibull or log-logistic.
#'
#' @param df Data frame with censoring columns.
#' @param distributions Character vector of `survreg` distributions to compare.
#'   Defaults to `c("lognormal", "weibull", "loglogistic")`.
#' @param col_detected Column name for detected values.
#' @param col_left_flag Logical column marking left-censored values.
#' @param col_left_limit Numeric LOD for left-censored values.
#' @param col_right_flag Logical column marking right-censored values.
#' @param col_right_limit Numeric upper limit for right-censored values.
#' @param min_total Minimum number of usable observations required to fit.
#' @param min_exact Minimum number of exact detections required to fit.
#'
#' @return A tibble with one row per distribution and columns for sample counts,
#'   log-likelihood, degrees of freedom, AIC, delta AIC, Akaike weight, best-fit
#'   flag, convergence flag, and any warning/error message.
#' @details
#' AIC comparison is not a formal proof that a distribution is true. It is a
#' likelihood-based diagnostic for whether the lognormal model fits better or
#' worse than alternative parametric censored models on the same data.
#' @export
compare_censored_distributions <- function(
    df,
    distributions = c("lognormal", "weibull", "loglogistic"),
    col_detected    = "detected_val",
    col_left_flag   = "left_cen",
    col_left_limit  = "left_limit",
    col_right_flag  = "right_cen",
    col_right_limit = "right_limit",
    min_total = 6,
    min_exact = 2
) {
  intervals <- .censored_intervals(
    df = df,
    col_detected = col_detected,
    col_left_flag = col_left_flag,
    col_left_limit = col_left_limit,
    col_right_flag = col_right_flag,
    col_right_limit = col_right_limit,
    caller = "compare_censored_distributions"
  )

  insufficient <- intervals$n_total < min_total || intervals$n_exact < min_exact

  fit_one <- function(dist) {
    supported <- c(
      "extreme", "logistic", "gaussian", "weibull", "exponential",
      "rayleigh", "loggaussian", "lognormal", "loglogistic", "t"
    )
    if (!dist %in% supported) {
      return(tibble::tibble(
        distribution = dist,
        n_total = intervals$n_total,
        n_left = intervals$n_left,
        n_right = intervals$n_right,
        n_exact = intervals$n_exact,
        logLik = NA_real_,
        df = NA_real_,
        AIC = NA_real_,
        converged = FALSE,
        message = paste0("unsupported distribution: ", dist)
      ))
    }

    if (insufficient) {
      return(tibble::tibble(
        distribution = dist,
        n_total = intervals$n_total,
        n_left = intervals$n_left,
        n_right = intervals$n_right,
        n_exact = intervals$n_exact,
        logLik = NA_real_,
        df = NA_real_,
        AIC = NA_real_,
        converged = FALSE,
        message = "insufficient_data"
      ))
    }

    warn <- character()
    fit <- tryCatch(
      withCallingHandlers(
        survival::survreg(
          survival::Surv(
            intervals$lower[intervals$ok],
            intervals$upper[intervals$ok],
            type = "interval2"
          ) ~ 1,
          dist = dist
        ),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(tibble::tibble(
        distribution = dist,
        n_total = intervals$n_total,
        n_left = intervals$n_left,
        n_right = intervals$n_right,
        n_exact = intervals$n_exact,
        logLik = NA_real_,
        df = NA_real_,
        AIC = NA_real_,
        converged = FALSE,
        message = conditionMessage(fit)
      ))
    }

    ll <- stats::logLik(fit)

    tibble::tibble(
      distribution = dist,
      n_total = intervals$n_total,
      n_left = intervals$n_left,
      n_right = intervals$n_right,
      n_exact = intervals$n_exact,
      logLik = as.numeric(ll),
      df = as.numeric(attr(ll, "df")),
      AIC = stats::AIC(fit),
      converged = TRUE,
      message = if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_
    )
  }

  out <- dplyr::bind_rows(lapply(distributions, fit_one))

  finite_aic <- is.finite(out$AIC)
  if (any(finite_aic)) {
    delta <- out$AIC - min(out$AIC[finite_aic])
    weights <- exp(-0.5 * delta)
    weights[!finite_aic] <- NA_real_
    weights <- weights / sum(weights, na.rm = TRUE)

    out$delta_AIC <- delta
    out$akaike_weight <- weights
    out$best <- finite_aic & out$delta_AIC == 0
  } else {
    out$delta_AIC <- NA_real_
    out$akaike_weight <- NA_real_
    out$best <- FALSE
  }

  out %>%
    dplyr::arrange(.data$AIC)
}
