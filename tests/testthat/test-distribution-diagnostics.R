test_that("compare_censored_distributions returns AIC table", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  out <- compare_censored_distributions(df)

  expect_equal(out$n_total, rep(8, 3))
  expect_setequal(out$distribution, c("lognormal", "weibull", "loglogistic"))
  expect_true(all(is.finite(out$AIC[out$converged])))
  expect_equal(sum(out$best), 1)
  expect_equal(sum(out$akaike_weight, na.rm = TRUE), 1, tolerance = 1e-8)
})

test_that("compare_censored_distributions reports insufficient data", {
  df <- tibble::tibble(
    detected_val = c(10, NA),
    left_cen = c(FALSE, TRUE),
    left_limit = c(NA_real_, 5),
    right_cen = c(FALSE, FALSE),
    right_limit = c(NA_real_, NA_real_)
  )

  out <- compare_censored_distributions(df)

  expect_true(all(!out$converged))
  expect_true(all(is.na(out$AIC)))
  expect_true(all(out$message == "insufficient_data"))
})

test_that("compare_censored_distributions captures fit errors", {
  df <- tibble::tibble(
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )

  out <- compare_censored_distributions(df, distributions = c("lognormal", "not_a_dist"))

  expect_true(out$converged[out$distribution == "lognormal"])
  expect_false(out$converged[out$distribution == "not_a_dist"])
  expect_true(is.na(out$AIC[out$distribution == "not_a_dist"]))
  expect_match(out$message[out$distribution == "not_a_dist"], "not_a_dist")
})
