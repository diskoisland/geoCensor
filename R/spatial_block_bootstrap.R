#' Experimental spatial block bootstrap for a declustered geometric mean
#'
#' Resamples occupied nonoverlapping square blocks with replacement and refits
#' a weighted censored lognormal model. This is an experimental sensitivity
#' procedure, not a validated confidence-interval method for general NGDB designs.
#'
#' @param df One medium/domain with standard censoring fields and coordinates.
#' @param x_col,y_col Numeric projected coordinate column names.
#' @param cell_size Declustering cell width, in coordinate units.
#' @param block_size Resampling block width, in the same units. This is distinct
#'   from cell_size; assess it against spatial dependence and domain size.
#' @param n_boot Number of replicates, an integer at least 99. Default 199 is
#'   intended for pilot work; larger runs are needed to assess tail stability.
#' @param conf_level Interval probability, strictly between zero and one.
#' @param n_offsets Declustering grid origins per axis.
#' @param block_origin Length-two finite vector specifying the fixed block origin.
#' @param min_blocks Minimum occupied blocks required, at least two. Default 10
#'   is a computational guard, not evidence of adequate spatial replication.
#' @param min_success Minimum successful replicate fraction, in (0, 1].
#' @param seed Optional nonnegative integer seed. Seeded calls restore the caller's
#'   random-number state; NULL uses and advances the current random-number stream.
#' @return List with summary, replicates, and blocks tibbles. The summary contains
#'   the declustered point estimate, experimental percentile limits, sample/block
#'   counts, success counts, and status. Limits are missing for insufficient
#'   blocks, an invalid point fit, or excessive failed replicates. Returned
#'   intervals always have status experimental_interval, not validated coverage.
#' @details
#' Usable chemistry and finite coordinates define a common analysis subset.
#' Each replicate draws as many blocks as are occupied. Declustering weights
#' are recalculated over the unique selected observations, then multiplied by
#' each block's draw count. This preserves multiplicity rather than letting
#' declustering of duplicate coordinates cancel repeated block draws.
#' Limits are empirical quantiles of successful log geometric means (type 7).
#' Failures are recorded and never replaced by the point estimate.
#'
#' The procedure conditions on sampled blocks and cannot recover unsampled
#' terrain or guarantee correction of informative sampling. Fixed block
#' boundaries, cross-block dependence, sparse blocks, and censoring can affect
#' coverage. Validate for the intended sampling design and assess multiple
#' block sizes/origins before interpreting intervals. This is not an
#' implementation of the generalised spatial bootstrap of Lahiri and Zhu.
#' @references
#' Lahiri, S. N. and Zhu, J. (2006). Resampling methods for spatial regression
#' models under a class of stochastic designs. Annals of Statistics 34, 1774-1813.
#' \doi{10.1214/009053606000000551}
#' @export
spatial_block_bootstrap <- function(df, x_col, y_col, cell_size, block_size,
                                     n_boot = 199, conf_level = 0.95, n_offsets = 10,
                                     block_origin = c(0, 0), min_blocks = 10,
                                     min_success = 0.9, seed = NULL) {
  whole <- function(x, lower) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= lower && x == floor(x)
  if (!whole(n_boot, 99)) stop("n_boot must be a whole number at least 99")
  if (!whole(min_blocks, 2)) stop("min_blocks must be a whole number at least two")
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size <= 0) stop("block_size must be positive")
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) stop("conf_level must be in (0, 1)")
  if (!is.numeric(min_success) || length(min_success) != 1L ||
      !is.finite(min_success) || min_success <= 0 || min_success > 1) stop("min_success must be in (0, 1]")
  if (!is.numeric(block_origin) || length(block_origin) != 2L ||
      any(!is.finite(block_origin))) stop("block_origin must contain two finite numbers")
  if (!is.null(seed) && (!whole(seed, 0) || seed > .Machine$integer.max)) stop("Invalid seed")
  intervals <- .censored_intervals(df, "detected_val", "left_cen", "left_limit",
                                    "right_cen", "right_limit", "spatial_block_bootstrap")
  w <- cell_declustering_weights(df, x_col, y_col, cell_size, n_offsets)
  usable <- intervals$ok & w > 0
  x <- df[usable, , drop = FALSE]
  x$.spatial_boot_weight <- cell_declustering_weights(x, x_col, y_col, cell_size, n_offsets)
  point <- fit_cen_lognorm(x, weight_col = ".spatial_boot_weight")
  bx <- floor((x[[x_col]] - block_origin[1]) / block_size)
  by <- floor((x[[y_col]] - block_origin[2]) / block_size)
  if (any(!is.finite(bx)) || any(!is.finite(by)) || any(abs(bx) > 2^52) || any(abs(by) > 2^52)) {
    stop("block_size is too small relative to coordinates and origin")
  }
  key <- paste(bx, by, sep = ":")
  id <- match(key, unique(key))
  k <- length(unique(id))
  first <- !duplicated(id)
  blocks <- tibble::tibble(block_id = seq_len(k), x_index = bx[first], y_index = by[first],
                           n_observations = tabulate(id, nbins = k))
  replicates <- tibble::tibble(replicate = seq_len(n_boot), geomean = NA_real_,
                               n_selected_blocks = NA_integer_, method = "not_run")
  status <- if (!is.finite(point$geomean)) "point_fit_unavailable" else if (k < min_blocks) "insufficient_blocks" else "experimental_interval"
  if (status == "experimental_interval") {
    if (!is.null(seed)) {
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit({
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
      }, add = TRUE)
      set.seed(seed)
    }
    for (i in seq_len(n_boot)) {
      draws <- tabulate(sample.int(k, k, replace = TRUE), nbins = k)
      weights <- .block_replicate_weights(x, id, draws, x_col, y_col, cell_size, n_offsets)
      replicate_data <- x
      replicate_data$.spatial_boot_weight <- weights
      fit <- fit_cen_lognorm(replicate_data, weight_col = ".spatial_boot_weight")
      replicates$geomean[i] <- fit$geomean
      replicates$method[i] <- fit$method
      replicates$n_selected_blocks[i] <- sum(draws > 0)
    }
  }
  success <- is.finite(replicates$geomean) & replicates$geomean > 0
  n_success <- sum(success)
  limits <- c(NA_real_, NA_real_)
  if (status == "experimental_interval") {
    if (n_success < ceiling(n_boot * min_success)) status <- "insufficient_successful_replicates"
    else limits <- exp(stats::quantile(log(replicates$geomean[success]),
                                       c((1 - conf_level) / 2, (1 + conf_level) / 2), names = FALSE))
  }
  summary <- tibble::tibble(
    geomean = point$geomean, ci_low = limits[1], ci_high = limits[2],
    conf_level = conf_level, cell_size = cell_size, block_size = block_size,
    origin_x = block_origin[1], origin_y = block_origin[2], n_origins = n_offsets^2,
    n_input = nrow(df), n_total = nrow(x), n_excluded = sum(!usable),
    n_left = point$n_left, n_right = point$n_right,
    n_exact = point$n_total - point$n_left - point$n_right,
    n_blocks = k, n_boot = n_boot, n_success = n_success,
    n_failed = sum(replicates$method != "not_run") - n_success,
    status = status, method = "experimental_nonoverlapping_block_percentile"
  )
  list(summary = summary, replicates = replicates, blocks = blocks)
}

.block_replicate_weights <- function(df, id, draws, x_col, y_col, cell_size, n_offsets) {
  selected <- draws[id] > 0
  weights <- numeric(nrow(df))
  weights[selected] <- cell_declustering_weights(df[selected, , drop = FALSE],
                                                  x_col, y_col, cell_size, n_offsets) * draws[id[selected]]
  weights / sum(weights)
}
