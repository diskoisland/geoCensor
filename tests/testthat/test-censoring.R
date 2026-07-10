test_that("parse_usgs_q handles numeric, suffix, and qualifier censoring", {
  out <- parse_usgs_q(
    qvalue = c(" -2 ", "10.1111", "5", "7", "8"),
    qualifier = c(" l ", NA, "", "G", ">")
  )

  expect_true(out$left_cen[1])
  expect_equal(out$left_limit[1], 2)
  expect_true(out$right_cen[2])
  expect_equal(out$right_limit[2], 10)
  expect_equal(out$detected_val[3], 5)
  expect_true(out$right_cen[4])
  expect_equal(out$right_limit[4], 7)
  expect_true(out$right_cen[5])
  expect_equal(out$right_limit[5], 8)
})

test_that("parse_usgs_q rejects mismatched qualifier lengths", {
  expect_error(
    parse_usgs_q(c("1", "2"), qualifier = "L"),
    "same length"
  )
})

test_that("fit_cen_lognorm estimates finite values with enough interval data", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  fit <- fit_cen_lognorm(df)

  expect_equal(fit$n_total, 8)
  expect_equal(fit$n_left, 2)
  expect_equal(fit$n_right, 1)
  expect_true(is.finite(fit$geomean))
  expect_true(fit$ci_low < fit$geomean)
  expect_true(fit$ci_high > fit$geomean)
})

test_that("fit_lognorm_with_quants returns fitted quantiles in mg/kg", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  fit <- fit_lognorm_with_quants(df)

  expect_equal(fit$method, "MLE")
  expect_true(is.finite(fit$geomean))
  expect_true(is.finite(fit$median))
  expect_true(fit$p90 > fit$median)
})
