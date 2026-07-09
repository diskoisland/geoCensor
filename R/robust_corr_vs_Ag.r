#' Robust bootstrap Spearman correlation of metals with silver
#'
#' Convenience wrapper around [robust_cor_table()] for the common case of
#' correlating a suite of metals against silver (`Ag_mgkg`). Each coefficient is
#' a Spearman rho with an empirical 95% confidence interval obtained from
#' `n_boot` bootstrap resamples of the complete (non-missing) paired
#' observations. Bootstrap resampling reduces sensitivity to leverage points,
#' skewed distributions, and extreme values that are common in soil and sediment
#' geochemistry, yielding coefficients that are more reproducible and less
#' influenced by outliers than conventional single-pass correlations. This is the
#' function used to build the element--silver correlation tables (e.g. Table S6
#' of the U.S. weather-modification silver inventory).
#'
#' @param df Wide metals data frame containing an `Ag_mgkg` column and one column
#'   per metal (e.g. from [build_metals_wide()] joined to the main table).
#' @param metals Character vector of metal column names to correlate with silver.
#' @param target Name of the silver (target) column (default `"Ag_mgkg"`).
#' @param n_boot Number of bootstrap resamples (default 1000).
#' @param min_pairs Minimum number of complete pairs required to compute rho
#'   (default 10).
#'
#' @return Tibble with columns `Element`, `N_pairs`, `Rho`, `CI_low_95`, and
#'   `CI_high_95`, sorted by descending rho.
#' @seealso [robust_cor_table()]
#' @export
robust_corr_vs_Ag <- function(df,
                              metals,
                              target    = "Ag_mgkg",
                              n_boot    = 1000,
                              min_pairs = 10) {
  robust_cor_table(
    df        = df,
    target    = target,
    metals    = metals,
    n_boot    = n_boot,
    min_pairs = min_pairs
  )
}
