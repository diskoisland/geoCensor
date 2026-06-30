# -------- PCA wrapper ----------

#' Run log10 PCA on metals with coverage / variance filters
#'
#' @param df Data frame (e.g. Idaho sediments or soils).
#' @param metal_cols Metal columns to consider.
#' @param min_non_na Minimum non-NA per metal to include.
#' @param min_complete Minimum complete rows required overall.
#'
#' @return List with pca (prcomp object), scores, used_cols, keep_rows.
#' @export
run_pca_robust <- function(df,
                           metal_cols,
                           min_non_na   = 200,
                           min_complete = 50) {
  x <- df %>% dplyr::select(dplyr::all_of(metal_cols))
  
  non_na_counts <- colSums(!is.na(x))
  keep_cols <- names(non_na_counts)[non_na_counts >= min_non_na]
  x <- x %>% dplyr::select(dplyr::all_of(keep_cols))
  
  ok <- stats::complete.cases(x)
  x  <- x[ok, , drop = FALSE]
  
  if (nrow(x) < min_complete) {
    stop("run_pca_robust: Not enough complete rows: ", nrow(x))
  }
  
  eps   <- 1e-6
  x_log <- log10(x + eps)
  
  sds <- apply(x_log, 2, stats::sd, na.rm = TRUE)
  non_const <- names(sds)[sds > 0]
  x_log2    <- x_log[, non_const, drop = FALSE]
  
  if (ncol(x_log2) < 2) {
    stop("run_pca_robust: fewer than 2 non-constant metals.")
  }
  
  pca <- stats::prcomp(x_log2, center = TRUE, scale. = TRUE)
  
  scores <- as.data.frame(pca$x) %>%
    dplyr::mutate(row_id = which(ok))
  
  list(
    pca       = pca,
    scores    = scores,
    used_cols = non_const,
    keep_rows = ok
  )
}