method_fit_sample <- function() {
  tibble::tibble(detected_val = c(1:8, 0), left_cen = FALSE, left_limit = NA_real_,
    right_cen = FALSE, right_limit = NA_real_,
    amethod = c("ES_SQ", "ES_H2O_P", "UNKNOWN", NA, rep("AES_HF", 4), "ES_SQ"))
}

method_summary_names <- c("n_semiquantitative_input", "n_semiquantitative_used",
  "n_method_unspecified_used", "n_method_not_reviewed_used",
  "review_semiquantitative", "method_review_status")

test_that("direct fits retain method counts for eligible observations and weights", {
  x <- method_fit_sample()
  fit <- fit_cen_lognorm(x)
  expect_equal(unlist(fit[1:3]), c(n_total = 8, n_left = 0, n_right = 0))
  expect_equal(fit$n_semiquantitative_input, 3)
  expect_equal(fit$n_semiquantitative_used, 2)
  expect_equal(fit$n_method_unspecified_used, 2)
  expect_equal(fit$n_method_not_reviewed_used, 4)
  expect_true(fit$review_semiquantitative)
  expect_equal(fit$method_review_status, "review_semiquantitative")
  expect_equal(fit$geomean, exp(mean(log(1:8))), tolerance = 1e-6)
  x$w <- c(0, 0, rep(1, 7))
  weighted <- fit_cen_lognorm(x, weight_col = "w")
  expect_equal(weighted$n_total, 6)
  expect_equal(weighted$n_semiquantitative_input, 3)
  expect_equal(weighted$n_semiquantitative_used, 0)
  expect_true(is.na(weighted$review_semiquantitative))
  expect_equal(weighted$method_review_status, "method_information_unresolved")
  expect_equal(weighted$geomean, exp(mean(log(3:8))), tolerance = 1e-6)
  expect_true(all(is.na(weighted[c("ci_low", "ci_high")])))
  names(x)[1:5] <- c("v", "lf", "ll", "rf", "rl")
  expect_equal(fit_cen_lognorm(x, "v", "lf", "ll", "rf", "rl", "w"), weighted)
})

test_that("method concerns remain visible when fits are unavailable", {
  x <- method_fit_sample()[1:5, ]
  fit <- fit_cen_lognorm(x)
  expect_equal(fit$method, "insufficient_data")
  expect_true(is.na(fit$geomean))
  expect_equal(fit$n_semiquantitative_used, 2)
  x$detected_val <- 0
  empty <- fit_cen_lognorm(x)
  expect_equal(empty$n_semiquantitative_input, 2)
  expect_equal(empty$n_semiquantitative_used, 0)
  expect_true(is.na(empty$review_semiquantitative))
  expect_equal(empty$method_review_status, "no_usable_data")
  zero <- fit_cen_lognorm(x[FALSE, ])
  expect_equal(zero$n_total, 0)
  expect_equal(zero$method_review_status, "no_usable_data")
  testthat::local_mocked_bindings(survreg = function(...) warning("failed"), .package = "survival")
  failed <- fit_cen_lognorm(method_fit_sample())
  expect_equal(failed$method, "failed_fit")
  expect_equal(failed$n_semiquantitative_used, 2)
  expect_true(failed$review_semiquantitative)
})

test_that("quantile and state summaries keep their own method-review subsets", {
  x <- method_fit_sample()
  q <- fit_lognorm_with_quants(x, probs = c(.25, .75))
  fit <- fit_cen_lognorm(x)
  expect_equal(q[method_summary_names], fit[method_summary_names])
  expect_equal(q$geomean, fit$geomean/1000)
  expect_equal(fit_lognorm_with_quants(x[1:5, ])$n_semiquantitative_used, 2)
  y <- x; y$amethod <- NA_character_; y$state <- "UT"; y$country <- "US"
  x$state <- "ID"; x$country <- "US"
  foreign <- x; foreign$country <- "CA"; foreign$state <- "BC"
  states <- mle_by_state(dplyr::bind_rows(x, y, foreign))
  expect_setequal(states$state, c("ID", "UT"))
  expect_equal(states$n_semiquantitative_used[states$state == "ID"], 2)
  expect_equal(states$n_semiquantitative_used[states$state == "UT"], 0)
  expect_equal(states$n_method_unspecified_used[states$state == "UT"], 8)
  expect_equal(states$method_review_status[states$state == "UT"], "method_information_unresolved")
})
