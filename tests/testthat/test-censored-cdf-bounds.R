bounds_sample <- function() {
  tibble::tibble(detected_val = c(3, NA_real_, NA_real_),
    left_cen = c(FALSE, TRUE, FALSE), left_limit = c(NA_real_, 2, NA_real_),
    right_cen = c(FALSE, FALSE, TRUE), right_limit = c(NA_real_, NA_real_, 4),
    w = c(.2, .5, .3))
}

test_that("censoring bounds contain every admissible empirical CDF", {
  x <- bounds_sample(); t <- c(1, 2, 3, 4, 5, 10)
  out <- censored_cdf_bounds(x, t, "w")
  expect_equal(out$cdf$cdf_lower, c(0, .5, .7, .7, .7, .7))
  expect_equal(out$cdf$cdf_upper, c(.5, .5, .7, .7, 1, 1))
  for (a in c(.25, 1, 2)) for (b in c(4.01, 5, 100)) {
    actual <- vapply(t, function(v) sum(x$w[c(3, a, b) <= v]), 0.0)
    expect_true(all(actual >= out$cdf$cdf_lower - 1e-14))
    expect_true(all(actual <= out$cdf$cdf_upper + 1e-14))
  }
  expect_equal(out$support$n_exact, 1)
  expect_equal(out$support$weighted_exact_fraction, .2)
})

test_that("exact ties and censoring-limit endpoints have correct inclusion", {
  x <- bounds_sample()
  x$detected_val <- c(3, 3, 5); x$left_cen[] <- x$right_cen[] <- FALSE
  out <- censored_cdf_bounds(x, c(5, 2, 3, 3), "w")
  expect_equal(out$cdf$threshold, c(2, 3, 5))
  expect_equal(out$cdf$cdf_lower, c(0, .7, 1))
  expect_equal(out$cdf$cdf_upper, out$cdf$cdf_lower)
  expect_equal(out$cdf$censoring_gap, rep(0, 3))
  right <- bounds_sample()[3, ]
  out <- censored_cdf_bounds(right, c(3, 4, 4.01))
  expect_equal(out$cdf$cdf_lower, c(0, 0, 0))
  expect_equal(out$cdf$cdf_upper, c(0, 0, 1))
  left <- bounds_sample()[2, ]
  out <- censored_cdf_bounds(left, c(1, 2, 3))
  expect_equal(out$cdf$cdf_lower, c(0, 1, 1))
  expect_equal(out$cdf$cdf_upper, c(1, 1, 1))
})

test_that("weight scaling and replicated records do not narrow censoring bounds", {
  x <- bounds_sample(); t <- 1:10
  ref <- censored_cdf_bounds(x, t, "w")
  x$w <- x$w * 1e300
  expect_equal(censored_cdf_bounds(x, t, "w"), ref)
  repeated <- x[rep(1:3, each = 5), ]
  expect_equal(censored_cdf_bounds(repeated, t, "w")$cdf, ref$cdf)
  x$w[] <- 1e308
  expect_equal(censored_cdf_bounds(x, t, "w")$cdf,
               censored_cdf_bounds(x, t)$cdf)
})

test_that("unusable chemistry and zero weights are excluded explicitly", {
  x <- bounds_sample(); x$left_limit[2] <- NA_real_; x$w[3] <- 0
  out <- censored_cdf_bounds(x, c(2, 3, 4), "w")
  expect_equal(out$support$n_excluded, 2)
  expect_equal(out$support$n_total, 1)
  expect_equal(out$cdf$cdf_lower, c(0, 1, 1))
  x$w[] <- 0
  empty <- censored_cdf_bounds(x, 1:3, "w")
  expect_equal(empty$support$n_total, 0)
  expect_true(all(is.na(empty$cdf$cdf_lower)))
  expect_true(is.na(empty$support$weighted_exact_fraction))
  expect_error(plot_censored_cdf_bounds(empty), "No usable")
  expect_equal(censored_cdf_bounds(x[FALSE, ], 1:3)$support$n_total, 0)
})

test_that("bounds validate arguments and support renamed columns", {
  x <- bounds_sample()
  expect_error(censored_cdf_bounds(x, c(0, 1)), "positive")
  expect_error(censored_cdf_bounds(x, c(1, Inf)), "positive")
  expect_error(censored_cdf_bounds(x, numeric()), "positive")
  expect_error(censored_cdf_bounds(x, 1, "missing"), "weight column")
  x$w[1] <- -1
  expect_error(censored_cdf_bounds(x, 1, "w"), "nonnegative")
  x$w[1] <- NA_real_
  expect_error(censored_cdf_bounds(x, 1, "w"), "finite")
  x <- bounds_sample(); names(x)[1:5] <- c("v", "l", "ll", "r", "rr")
  out <- censored_cdf_bounds(x, 1:5, col_detected = "v", col_left_flag = "l",
    col_left_limit = "ll", col_right_flag = "r", col_right_limit = "rr")
  expect_equal(out$cdf, censored_cdf_bounds(bounds_sample(), 1:5)$cdf)
  p <- plot_censored_cdf_bounds(out, title = "Defined domain")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Defined domain")
  expect_match(p$labels$caption, "sampling and spatial")
})

test_that("CDF support retains method counts on its positive-weight valid subset", {
  x <- bounds_sample()
  x$amethod <- c(" ES_SQ ", "es_h2o_p", "UNKNOWN")
  x <- dplyr::bind_rows(x, tibble::tibble(detected_val = 0, left_cen = FALSE,
    left_limit = NA_real_, right_cen = FALSE, right_limit = NA_real_, w = 1,
    amethod = "ES_SQ"))
  out <- censored_cdf_bounds(x, 1:5, "w")
  expect_equal(out$support$n_semiquantitative_input, 3)
  expect_equal(out$support$n_semiquantitative_used, 2)
  expect_equal(out$support$n_method_unspecified_used, 1)
  expect_equal(out$support$n_method_not_reviewed_used, 0)
  expect_true(out$support$review_semiquantitative)
  expect_equal(out$support$method_review_status, "review_semiquantitative")
  expect_equal(out$cdf, censored_cdf_bounds(bounds_sample(), 1:5, "w")$cdf)
  x$w[1] <- 0; x$left_limit[2] <- NA_real_
  restricted <- censored_cdf_bounds(x, 1:5, "w")
  expect_equal(restricted$support$n_total, 1)
  expect_equal(restricted$support$n_semiquantitative_input, 3)
  expect_equal(restricted$support$n_semiquantitative_used, 0)
  expect_true(is.na(restricted$support$review_semiquantitative))
  expect_equal(restricted$support$method_review_status, "method_information_unresolved")
  names(x)[1:5] <- c("v", "lf", "ll", "rf", "rl")
  renamed <- censored_cdf_bounds(x, 1:5, "w", "v", "lf", "ll", "rf", "rl")
  expect_equal(renamed, restricted)
})

test_that("coincident CDF bounds do not clear semiquantitative measurement concerns", {
  x <- bounds_sample()
  x$detected_val <- c(3, 3, 5); x$left_cen[] <- x$right_cen[] <- FALSE
  x$amethod <- "ES_SQ"
  out <- censored_cdf_bounds(x, 1:5)
  expect_true(all(out$cdf$censoring_gap == 0))
  expect_true(out$support$review_semiquantitative)
  expect_equal(out$support$n_semiquantitative_used, 3)
  p <- plot_censored_cdf_bounds(out)
  expect_match(p$labels$caption, "Known semiquantitative codes: 3/3 usable records", fixed = TRUE)
  expect_match(p$labels$caption, "measurement precision unresolved", fixed = TRUE)
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("CDF method information remains unresolved for missing and empty data", {
  x <- bounds_sample()
  absent <- censored_cdf_bounds(x, 1:5)
  expect_equal(absent$support$n_method_unspecified_used, 3)
  expect_true(is.na(absent$support$review_semiquantitative))
  expect_equal(absent$support$method_review_status, "method_information_unresolved")
  x$amethod <- "ES_SQ"; x$w <- 0
  empty <- censored_cdf_bounds(x, 1:5, "w")
  expect_equal(empty$support$n_semiquantitative_input, 3)
  expect_equal(empty$support$n_semiquantitative_used, 0)
  expect_equal(empty$support$method_review_status, "no_usable_data")
  expect_true(is.na(empty$support$review_semiquantitative))
  expect_true(all(is.na(empty$cdf$cdf_lower)))
  expect_equal(censored_cdf_bounds(x[FALSE, ], 1:5)$support$method_review_status, "no_usable_data")
})

test_that("older saved CDF results plot without invented method counts", {
  out <- censored_cdf_bounds(bounds_sample(), 1:5)
  method_fields <- c("n_semiquantitative_input", "n_semiquantitative_used",
    "n_method_unspecified_used", "n_method_not_reviewed_used",
    "review_semiquantitative", "method_review_status")
  legacy <- out
  legacy$support <- legacy$support[, setdiff(names(legacy$support), method_fields)]
  legacy <- unserialize(serialize(legacy, NULL))
  p <- plot_censored_cdf_bounds(legacy)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "Analytical-method review unavailable in this saved result", fixed = TRUE)
  expect_false(grepl("Known semiquantitative codes:", p$labels$caption, fixed = TRUE))
  expect_silent(ggplot2::ggplot_build(p))
})
