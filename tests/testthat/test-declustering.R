clustered_sample <- function() {
  tibble::tibble(
    x = c(rep(0, 8), 100, 200), y = 0,
    detected_val = c(seq(80, 150, by = 10), 10, 20),
    left_cen = FALSE, left_limit = NA_real_,
    right_cen = FALSE, right_limit = NA_real_
  )
}

test_that("occupied cells share influence across repeated locations", {
  x <- clustered_sample()
  w <- cell_declustering_weights(x, "x", "y", 20)
  expect_equal(sum(w), 1)
  expect_equal(sum(w[1:8]), 1 / 3)
  expect_equal(w[9:10], c(1 / 3, 1 / 3))
  expect_equal(w, cell_declustering_weights(x, "x", "y", 20))
  shifted <- dplyr::mutate(x, x = .data$x + 500, y = .data$y - 200)
  expect_equal(w, cell_declustering_weights(shifted, "x", "y", 20))
  scaled <- dplyr::mutate(x, x = .data$x / 1000, y = .data$y / 1000)
  expect_equal(w, cell_declustering_weights(scaled, "x", "y", 0.02))
})

test_that("invalid coordinates get zero weight and empty inputs are supported", {
  x <- clustered_sample()
  x$x[1] <- NA_real_
  expect_equal(cell_declustering_weights(x, "x", "y", 20)[1], 0)
  expect_length(cell_declustering_weights(x[FALSE, ], "x", "y", 20), 0)
  x$x[] <- NA_real_
  expect_equal(cell_declustering_weights(x, "x", "y", 20), rep(0, 10))
  expect_error(cell_declustering_weights(x, "x", "y", 0), "positive")
  expect_error(cell_declustering_weights(x, "x", "y", 20, 1.5), "whole number")
})

test_that("weighted estimates match weighted log means without ordinary intervals", {
  x <- clustered_sample()
  x$w <- cell_declustering_weights(x, "x", "y", 20)
  fit <- fit_cen_lognorm(x, weight_col = "w")
  expect_equal(fit$geomean, exp(stats::weighted.mean(log(x$detected_val), x$w)), tolerance = 1e-6)
  expect_lt(fit$geomean, fit_cen_lognorm(x)$geomean)
  expect_true(is.na(fit$ci_low) && is.na(fit$ci_high))
  expect_equal(fit$method, "weighted_survreg_lognormal")
  x$w <- x$w * 1000
  expect_equal(fit, fit_cen_lognorm(x, weight_col = "w"))
  x$w[1] <- -1
  expect_error(fit_cen_lognorm(x, weight_col = "w"), "nonnegative")
})

test_that("weighted likelihood preserves both censoring directions", {
  x <- clustered_sample()
  x$left_cen[9] <- TRUE
  x$left_limit[9] <- 10
  x$right_cen[8] <- TRUE
  x$right_limit[8] <- 150
  x$w <- cell_declustering_weights(x, "x", "y", 20)
  fit <- fit_cen_lognorm(x, weight_col = "w")
  expect_equal(fit$n_left, 1)
  expect_equal(fit$n_right, 1)
  lower <- ifelse(x$left_cen, 0, ifelse(x$right_cen, x$right_limit, x$detected_val))
  upper <- ifelse(x$left_cen, x$left_limit, ifelse(x$right_cen, Inf, x$detected_val))
  reference <- survival::survreg(survival::Surv(lower, upper, type = "interval2") ~ 1,
                                weights = x$w * nrow(x), dist = "lognormal")
  expect_equal(fit$geomean, unname(exp(stats::coef(reference)[1])), tolerance = 1e-6)
  x$w[9] <- 0
  expect_equal(fit_cen_lognorm(x, weight_col = "w")$n_left, 0)
})

test_that("sensitivity uses a common valid subset and keeps all requested scales", {
  x <- clustered_sample()
  x$x[1] <- NA_real_
  x$detected_val[2] <- NA_real_
  out <- declustering_sensitivity(x, c(10, 20, 40), "x", "y")
  expect_equal(out$cell_size, c(10, 20, 40))
  expect_equal(out$n_total, rep(8, 3))
  expect_equal(out$n_excluded, rep(2, 3))
  expect_equal(out$n_origins, rep(100, 3))
  expect_equal(out$exact_fraction, rep(1, 3))
  expect_equal(out$weighted_exact_fraction, rep(1, 3))
  expect_equal(out$weighted_left_fraction, rep(0, 3))
  expect_equal(out$weighted_right_fraction, rep(0, 3))
  expect_true(all(is.finite(out$geomean)))
  expect_true(all(out$weight_effective_n < out$n_total))
  expect_equal(out$unweighted_geomean, rep(fit_cen_lognorm(x[3:10, ])$geomean, 3))
  x$x[] <- NA_real_
  empty <- declustering_sensitivity(x, 20, "x", "y")
  expect_equal(empty$n_total, 0)
  expect_equal(empty$method, "insufficient_data")
  expect_true(all(is.na(empty[, c("exact_fraction", "weighted_exact_fraction",
    "weighted_left_fraction", "weighted_right_fraction")])))
  expect_error(declustering_sensitivity(x, numeric(), "x", "y"), "positive")
})

test_that("declustering reports censoring weight rather than raw record fractions", {
  x <- clustered_sample()
  x$left_cen[9] <- TRUE
  x$left_limit[9] <- 10
  x$right_cen[10] <- TRUE
  x$right_limit[10] <- 20
  out <- declustering_sensitivity(x, c(10, 20), "x", "y")
  # Eight exact values occupy one cell; each isolated censoring type another.
  expect_equal(out$exact_fraction, rep(.8, 2))
  expect_equal(out$weighted_exact_fraction, rep(1 / 3, 2))
  expect_equal(out$weighted_left_fraction, rep(1 / 3, 2))
  expect_equal(out$weighted_right_fraction, rep(1 / 3, 2))
  expect_equal(out$weighted_exact_fraction + out$weighted_left_fraction +
                 out$weighted_right_fraction, rep(1, 2))

  x$x[10] <- NA_real_
  x$detected_val[1] <- NA_real_
  excluded <- declustering_sensitivity(x, 20, "x", "y")
  expect_equal(excluded$n_excluded, 2)
  expect_equal(excluded$exact_fraction, 7 / 8)
  expect_equal(excluded$weighted_exact_fraction, .5)
  expect_equal(excluded$weighted_left_fraction, .5)
  expect_equal(excluded$weighted_right_fraction, 0)

  x <- clustered_sample()
  x$left_cen[] <- TRUE
  x$left_limit[] <- 20
  censored <- declustering_sensitivity(x, 20, "x", "y")
  expect_equal(censored$method, "insufficient_data")
  expect_equal(censored$weighted_exact_fraction, 0)
  expect_equal(censored$weighted_left_fraction, 1)
  expect_equal(censored$weighted_right_fraction, 0)
})
