test_that("plot_lognormal_fit returns a ggplot object", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  p <- plot_lognormal_fit(df, title = "Test fit")

  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 5)
  expect_equal(p$labels$title, "Test fit")
  expect_match(p$labels$caption, "GM = exp\\(meanlog\\)")
  expect_match(p$labels$caption, "exact detections only")
})

test_that("plot_lognormal_fit can hide detected points in CDF view", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  p <- plot_lognormal_fit(df, show_detects = FALSE)

  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 4)
})

test_that("plot_lognormal_fit can show nonparametric censored CDF", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  p <- plot_lognormal_fit(df, show_nonparametric = TRUE)

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "nonparametric censored CDF")
})

test_that("plot_lognormal_fit still supports density view", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  p <- plot_lognormal_fit(df, plot_type = "density")

  expect_s3_class(p, "ggplot")
  expect_equal(
    p$labels$caption,
    "Histogram shows exact detections only; fitted curve uses exact and censored observations."
  )
})

test_that("plot_lognormal_fit supports custom censoring column names", {
  df <- tibble::tibble(
    val = c(NA, NA, 1, 2, 3, 4, 5, 6),
    lflag = c(TRUE, TRUE, rep(FALSE, 6)),
    llim = c(0.4, 0.6, rep(NA_real_, 6)),
    rflag = rep(FALSE, 8),
    rlim = rep(NA_real_, 8)
  )

  p <- plot_lognormal_fit(
    df,
    col_detected = "val",
    col_left_flag = "lflag",
    col_left_limit = "llim",
    col_right_flag = "rflag",
    col_right_limit = "rlim"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_lognormal_fit validates minimum data requirements", {
  df <- tibble::tibble(
    detected_val = c(1, NA),
    left_cen = c(FALSE, TRUE),
    left_limit = c(NA_real_, 0.5),
    right_cen = c(FALSE, FALSE),
    right_limit = c(NA_real_, NA_real_)
  )

  expect_error(
    plot_lognormal_fit(df),
    "insufficient data"
  )
})
