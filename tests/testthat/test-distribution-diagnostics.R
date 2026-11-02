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

distribution_method_sample <- function() {
  tibble::tibble(detected_val = c(20,35,60,80,110,200,300,500,0), left_cen = FALSE,
    left_limit = NA_real_, right_cen = FALSE, right_limit = NA_real_,
    amethod = c("ES_SQ", "ES_H2O_P", NA, "UNKNOWN", rep("AES_HF", 4), "ES_SQ"))
}

test_that("distribution comparisons retain method review independently of AIC", {
  x <- distribution_method_sample()
  out <- compare_censored_distributions(x)
  expect_equal(out$n_total, rep(8, 3))
  expect_equal(out$n_semiquantitative_input, rep(3, 3))
  expect_equal(out$n_semiquantitative_used, rep(2, 3))
  expect_equal(out$n_method_unspecified_used, rep(2, 3))
  expect_equal(out$n_method_not_reviewed_used, rep(4, 3))
  expect_true(all(out$review_semiquantitative))
  expect_true(out$review_semiquantitative[out$best])
  expect_true(all(out$method_review_status == "review_semiquantitative"))
  expect_equal(anyDuplicated(names(out)), 0L)
  reference <- compare_censored_distributions(x[, setdiff(names(x), "amethod")])
  numeric_columns <- c("distribution", "logLik", "df", "AIC", "converged", "delta_AIC", "akaike_weight", "best")
  expect_equal(out[numeric_columns], reference[numeric_columns])
  names(x)[1:5] <- c("v", "lf", "ll", "rf", "rl")
  expect_equal(compare_censored_distributions(x, col_detected = "v", col_left_flag = "lf",
    col_left_limit = "ll", col_right_flag = "rf", col_right_limit = "rl"), out)
})

test_that("insufficient, unsupported and failed distribution fits retain method concerns", {
  x <- distribution_method_sample()
  sparse <- compare_censored_distributions(x[1:5, ])
  expect_true(all(sparse$message == "insufficient_data"))
  expect_true(all(is.na(sparse$AIC)))
  expect_equal(sparse$n_semiquantitative_used, rep(2, 3))
  unsupported <- compare_censored_distributions(x, distributions = "unsupported")
  expect_true(unsupported$review_semiquantitative)
  expect_false(unsupported$best)
  original <- survival::survreg
  testthat::local_mocked_bindings(survreg = function(...) {
    warning("did not converge"); original(...)
  }, .package = "survival")
  failed <- compare_censored_distributions(x)
  expect_true(all(failed$review_semiquantitative))
  expect_true(all(!failed$converged))
  expect_true(all(!failed$best))
  expect_true(all(is.na(failed$akaike_weight)))
  expect_equal(failed$n_semiquantitative_used, rep(2, 3))
})

test_that("empty or unknown method information never creates a clear distribution flag", {
  x <- distribution_method_sample()
  x$amethod <- NA_character_
  absent <- compare_censored_distributions(x)
  expect_true(all(is.na(absent$review_semiquantitative)))
  expect_equal(absent$n_method_unspecified_used, rep(8, 3))
  expect_true(all(absent$method_review_status == "method_information_unresolved"))
  x$detected_val <- 0; x$amethod <- "ES_SQ"
  invalid <- compare_censored_distributions(x)
  expect_true(all(invalid$method_review_status == "no_usable_data"))
  expect_equal(invalid$n_semiquantitative_input, rep(9, 3))
  expect_equal(invalid$n_semiquantitative_used, rep(0, 3))
  expect_true(all(is.na(invalid$review_semiquantitative)))
  empty <- compare_censored_distributions(x[FALSE, ])
  expect_true(all(empty$method_review_status == "no_usable_data"))
})
