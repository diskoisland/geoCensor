bootstrap_sample <- function() {
  grid <- expand.grid(x = 0:3, y = 0:3)
  x <- grid[rep(seq_len(nrow(grid)), each = 3), ]
  x$x <- x$x + 0.2
  x$y <- x$y + 0.2
  x$detected_val <- exp(rep(seq(2, 4, length.out = 16), each = 3) + rep(c(-0.1, 0, 0.1), 16))
  x$left_cen <- FALSE; x$left_limit <- NA_real_
  x$right_cen <- FALSE; x$right_limit <- NA_real_
  x
}

test_that("block draws preserve multiplicity after weights are recalculated", {
  x <- bootstrap_sample()
  id <- rep(1:16, each = 3)
  draws <- c(2L, 0L, rep(1L, 14))
  w <- .block_replicate_weights(x, id, draws, "x", "y", 1, 2)
  expect_equal(sum(w[id == 1]), 2 / 16)
  expect_equal(sum(w[id == 2]), 0)
  expect_equal(sum(w), 1)
})

test_that("bootstrap is reproducible, preserves seed, and returns empirical limits", {
  x <- bootstrap_sample()
  x$left_cen[1:3] <- TRUE; x$left_limit[1:3] <- 10
  x$right_cen[46:48] <- TRUE; x$right_limit[46:48] <- 50
  set.seed(42); previous <- .Random.seed
  a <- spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 99, n_offsets = 2, seed = 8)
  expect_identical(.Random.seed, previous)
  b <- spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 99, n_offsets = 2, seed = 8)
  expect_equal(a, b)
  expect_equal(a$summary$status, "experimental_interval")
  expect_equal(a$summary$n_blocks, 16)
  expect_equal(a$summary$n_left, 3)
  expect_equal(a$summary$n_right, 3)
  success <- a$replicates$geomean[is.finite(a$replicates$geomean)]
  expect_equal(c(a$summary$ci_low, a$summary$ci_high),
               exp(stats::quantile(log(success), c(0.025, 0.975), names = FALSE)))
  expect_equal(a$summary$n_success + a$summary$n_failed, 99)
})

test_that("too few blocks and unusable inputs suppress intervals", {
  x <- bootstrap_sample()
  a <- spatial_block_bootstrap(x, "x", "y", 1, 10, n_boot = 99, n_offsets = 2)
  expect_equal(a$summary$status, "insufficient_blocks")
  expect_true(is.na(a$summary$ci_low))
  expect_equal(a$summary$n_success, 0)
  expect_true(all(a$replicates$method == "not_run"))
  x$x[] <- NA_real_
  empty <- spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 99, n_offsets = 2)
  expect_equal(empty$summary$n_total, 0)
  expect_equal(empty$summary$status, "point_fit_unavailable")
})

test_that("failed bootstrap fits remain visible and disable intervals", {
  x <- bootstrap_sample()
  original <- fit_cen_lognorm
  calls <- 0
  testthat::local_mocked_bindings(fit_cen_lognorm = function(...) {
    calls <<- calls + 1
    out <- original(...)
    if (calls > 1) {out$geomean <- NA_real_; out$method <- "failed_fit"}
    out
  })
  a <- spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 99, n_offsets = 1, seed = 1)
  expect_equal(a$summary$status, "insufficient_successful_replicates")
  expect_equal(a$summary$n_failed, 99)
  expect_true(is.na(a$summary$ci_low))
  expect_true(all(a$replicates$method == "failed_fit"))
})

test_that("bootstrap validates controls", {
  x <- bootstrap_sample()
  expect_error(spatial_block_bootstrap(x, "x", "y", 1, 0), "block_size")
  expect_error(spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 20), "99")
  expect_error(spatial_block_bootstrap(x, "x", "y", 1, 1, conf_level = 1), "conf_level")
  expect_error(spatial_block_bootstrap(x, "x", "y", 1, 1, block_origin = NA_real_), "block_origin")
  expect_error(spatial_block_bootstrap(x, "x", "y", 1, 1, seed = -1), "seed")
})

test_that("uncensored block bootstrap equals resampling site log means", {
  x <- bootstrap_sample()
  id <- rep(1:16, each = 3)
  site_means <- as.numeric(tapply(log(x$detected_val), id, mean))
  set.seed(193)
  expected <- replicate(99, mean(site_means[sample.int(16, 16, replace = TRUE)]))
  out <- spatial_block_bootstrap(x, "x", "y", 1, 1, n_boot = 99,
                                  n_offsets = 2, seed = 193)
  expect_equal(log(out$summary$geomean), mean(site_means), tolerance = 1e-7)
  expect_equal(log(out$replicates$geomean), expected, tolerance = 1e-7)
  # The conditional bootstrap variance of a sample mean is known exactly.
  draws <- lapply(seq_len(16), function(i) {z <- integer(16); z[i] <- 16L; z})
  for (i in c(1, 8, 16)) {
    w <- .block_replicate_weights(x, id, draws[[i]], "x", "y", 1, 2)
    expect_equal(sum(w[id == i]), 1)
    expect_equal(sum(w[id != i]), 0)
  }
})
