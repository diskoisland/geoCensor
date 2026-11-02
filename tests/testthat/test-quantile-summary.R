quantile_sample <- function() {
  tibble::tibble(
    detected_val = c(20, 35, 60, 80, 110, 200, 300, 500),
    left_cen = FALSE, left_limit = NA_real_,
    right_cen = FALSE, right_limit = NA_real_
  )
}

test_that("uncensored quantile summaries match the analytic lognormal MLE", {
  x <- quantile_sample()
  mu <- mean(log(x$detected_val))
  sigma <- sqrt(mean((log(x$detected_val) - mu)^2))
  fit <- fit_lognorm_with_quants(x)
  expect_equal(fit$method, "MLE")
  expect_equal(fit$n_left, 0)
  expect_equal(fit$meanlog, mu, tolerance = 1e-6)
  expect_equal(fit$sdlog, sigma, tolerance = 1e-6)
  expect_equal(fit$geomean, exp(mu)/1000, tolerance = 1e-6)
  expect_equal(fit$p90, stats::qlnorm(.9, mu, sigma)/1000, tolerance = 1e-6)
  expect_equal(fit$median, fit$geomean)
  expect_named(fit, c("method", "n_total", "n_left", "n_right", "meanlog", "sdlog",
    "geomean", "ci_low", "ci_high", "p10", "median", "p90", "p95", "p99"))
})

test_that("right-only censoring and quantile units agree with the main fit", {
  x <- quantile_sample()
  x$right_cen[8] <- TRUE
  x$right_limit[8] <- 500
  fit <- fit_lognorm_with_quants(x)
  main <- fit_cen_lognorm(x)
  expect_equal(fit$n_right, 1)
  expect_equal(fit$geomean, main$geomean/1000)
  expect_equal(fit$ci_low, main$ci_low/1000)
  expect_equal(fit$ci_high, main$ci_high/1000)
})

test_that("quantile summaries exclude invalid limits and reject contradictory flags", {
  x <- quantile_sample()
  x$left_cen[1] <- TRUE
  fit <- fit_lognorm_with_quants(x)
  expect_equal(fit$n_total, 7)
  expect_equal(fit$geomean, exp(mean(log(x$detected_val[-1])))/1000, tolerance = 1e-6)
  x$right_cen[1] <- TRUE
  expect_error(fit_lognorm_with_quants(x), "both left- and right-censored")
  expect_error(fit_lognorm_with_quants(x[, -1]), "Missing columns")
})

test_that("insufficient and failed quantile summaries leave every requested value missing", {
  sparse <- fit_lognorm_with_quants(quantile_sample()[1:5, ], probs = c(.25, .75))
  expect_equal(sparse$method, "insufficient_data")
  expect_true(all(is.na(sparse[c("meanlog", "sdlog", "geomean", "ci_low", "ci_high", "median", "p25", "p75")])))
  one_exact <- quantile_sample()
  one_exact$left_cen[1:7] <- TRUE
  one_exact$left_limit[1:7] <- 50
  expect_equal(fit_lognorm_with_quants(one_exact)$method, "insufficient_data")
  expect_true(is.na(fit_lognorm_with_quants(one_exact)$median))
  original <- survival::survreg
  testthat::local_mocked_bindings(survreg = function(...) {
    warning("did not converge")
    original(...)
  }, .package = "survival")
  failed <- fit_lognorm_with_quants(quantile_sample(), probs = c(.25, .75))
  expect_equal(failed$method, "failed_fit")
  expect_true(all(is.na(failed[c("meanlog", "sdlog", "geomean", "ci_low", "ci_high", "median", "p25", "p75")])))
})

test_that("custom probabilities are returned without losing the model median", {
  x <- quantile_sample()
  mu <- mean(log(x$detected_val))
  sigma <- sqrt(mean((log(x$detected_val) - mu)^2))
  fit <- fit_lognorm_with_quants(x, probs = c(.25, .75))
  expect_equal(fit$p25, stats::qlnorm(.25, mu, sigma)/1000, tolerance = 1e-6)
  expect_equal(fit$p75, stats::qlnorm(.75, mu, sigma)/1000, tolerance = 1e-6)
  expect_equal(fit$median, exp(mu)/1000, tolerance = 1e-6)
  endpoints <- fit_lognorm_with_quants(x, probs = c(0, 1))
  expect_equal(endpoints$p0, 0)
  expect_equal(endpoints$p100, Inf)
  for (bad in list(numeric(), c(.5, .5), NA_real_, NaN, Inf, -0.1, 1.1, "0.5")) {
    expect_error(fit_lognorm_with_quants(x, probs = bad), "unique finite numeric probabilities")
  }
})
