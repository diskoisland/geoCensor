#' Robust Spearman correlation (bootstrap CI) of metals vs a target
#'
#' @param df Wide metals data (e.g. from build_metals_wide + join_main).
#' @param target Name of target column (e.g. "Ag_mgkg").
#' @param metals Vector of metal column names.
#' @param n_boot Number of bootstrap resamples.
#' @param min_pairs Minimum non-NA pairs required to compute rho.
#'
#' @return Tibble with Element, N_pairs, Rho, CI_low_95, CI_high_95.
#' @export
robust_cor_table <- function(df,
                             target,
                             metals,
                             n_boot   = 1000,
                             min_pairs = 50) {
  out_list <- lapply(metals, function(m) {
    x <- df[[target]]
    y <- df[[m]]

    ok <- is.finite(x) & is.finite(y)
    n  <- sum(ok)

    if (n < min_pairs) {
      return(tibble::tibble(
        Element    = m,
        N_pairs    = n,
        Rho        = NA_real_,
        CI_low_95  = NA_real_,
        CI_high_95 = NA_real_
      ))
    }

    rho <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))

    if (!is.finite(rho)) {
      return(tibble::tibble(
        Element    = m,
        N_pairs    = n,
        Rho        = NA_real_,
        CI_low_95  = NA_real_,
        CI_high_95 = NA_real_
      ))
    }

    boot_rho <- replicate(n_boot, {
      idx <- sample(which(ok), size = n, replace = TRUE)
      stats::cor(x[idx], y[idx], method = "spearman")
    })

    tibble::tibble(
      Element    = m,
      N_pairs    = n,
      Rho        = rho,
      CI_low_95  = stats::quantile(boot_rho, 0.025, na.rm = TRUE),
      CI_high_95 = stats::quantile(boot_rho, 0.975, na.rm = TRUE)
    )
  })

  dplyr::bind_rows(out_list) %>%
    dplyr::arrange(dplyr::desc(.data$Rho))
}
