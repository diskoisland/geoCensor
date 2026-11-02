exact_sample <- function() {
  tibble::tibble(
    detected_val = c(20, 35, 60, 80, 110, 200, 300, 500),
    left_cen = FALSE, left_limit = NA_real_,
    right_cen = FALSE, right_limit = NA_real_
  )
}

test_that("uncensored and right-censored samples can be estimated", {
  x <- exact_sample()
  fit <- fit_cen_lognorm(x)
  expect_equal(fit$geomean, exp(mean(log(x$detected_val))), tolerance = 1e-6)
  x$right_cen[8] <- TRUE
  x$right_limit[8] <- 500
  fit <- fit_cen_lognorm(x)
  expect_equal(fit$n_right, 1)
  expect_true(is.finite(fit$geomean))
})

test_that("invalid censoring limits never turn censored rows into detections", {
  x <- exact_sample()
  x$left_cen[1] <- TRUE
  x$left_limit[1] <- -1
  expect_equal(fit_cen_lognorm(x)$n_total, 7)
  expect_equal(unique(compare_censored_distributions(x)$n_exact), 7)
  intervals <- plot_censored_intervals(
    x, "detected_val", "left_cen", "left_limit", "right_cen", "right_limit", "test"
  )
  expect_false(intervals$ok[1])
  x$right_cen[1] <- TRUE
  expect_error(fit_cen_lognorm(x), "both left- and right-censored")
  expect_error(compare_censored_distributions(x), "both left- and right-censored")
})

test_that("failed convergence is excluded from estimates and AIC ranking", {
  original <- survival::survreg
  testthat::local_mocked_bindings(
    survreg = function(...) {
      warning("Ran out of iterations and did not converge")
      original(...)
    },
    .package = "survival"
  )
  x <- exact_sample()
  expect_equal(fit_cen_lognorm(x)$method, "failed_fit")
  expect_true(is.na(fit_cen_lognorm(x)$geomean))
  out <- compare_censored_distributions(x)
  expect_true(all(!out$converged))
  expect_true(all(!out$best))
  expect_true(all(is.na(out$akaike_weight)))
  expect_true(all(grepl("did not converge", out$message)))
})

test_that("distribution comparisons reject ambiguous arguments", {
  x <- exact_sample()
  expect_error(compare_censored_distributions(x, character()), "nonempty")
  expect_error(compare_censored_distributions(x, c("weibull", "weibull")), "unique")
  expect_error(compare_censored_distributions(x, NA_character_), "unique")
  expect_error(compare_censored_distributions(x, min_total = NA_real_), "whole numbers")
  expect_error(compare_censored_distributions(x, min_exact = 1.5), "whole numbers")
})
