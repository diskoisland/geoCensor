upper_recovery_input <- function(method = "CB_IRC", summary = "10.11111, CB_IRC",
                                 value = "10.11", qualifier = "") {
  tibble::tibble(lab_id = "001", species = "S", unit = "pct", qvalue = value,
    amethod = method, summary = summary, qualifier = qualifier)
}

test_that("summary recovery requires a known matching method and retains originals", {
  x <- dplyr::bind_rows(upper_recovery_input(), upper_recovery_input(" cb_irc "),
    upper_recovery_input("AES_HF"), upper_recovery_input(NA_character_),
    upper_recovery_input(""), upper_recovery_input("UNKNOWN"),
    upper_recovery_input("CB_IRC", "10.11111, UNKNOWN"))
  x$lab_id <- sprintf("%03d", seq_len(nrow(x)))
  out <- build_element(x, "S", to_ugkg = FALSE)
  expect_identical(out[names(x)], x)
  expect_identical(out$qvalue_raw, x$qvalue)
  expect_equal(out$right_cen, c(TRUE, TRUE, rep(FALSE, 5)))
  expect_equal(out$right_limit[1:2], c(10, 10))
  expect_equal(out$detected_val[3:7], rep(10.11, 5))
  expect_equal(out$censoring_source, c(rep("summary_upper_limit", 2), rep("qvalue", 5)))
  expect_equal(out$summary_recovery_status, c(rep("recovered_summary_upper_limit", 2),
    "blocked_method_mismatch", rep("blocked_method_information_missing", 4)))
  missing_column <- x[1, ]; missing_column$amethod <- NULL
  expect_equal(build_element(missing_column, "S")$summary_recovery_status, "blocked_method_information_missing")
})

test_that("only a complete first entry can supply an upper limit", {
  summaries <- c("10.11111, CB_IRC; 20, AES_HF", "20, AES_HF; 10.11111, CB_IRC",
    "10.11111, AES_HF; 10.11111, CB_IRC", "10.11111, CB_IRC, junk",
    "10.11111, CB_IRC,", "10.11111", "10.11111,", NA_character_, "", "bad, CB_IRC")
  x <- upper_recovery_input()[rep(1, length(summaries)), ]; x$summary <- summaries
  out <- build_element(x, "S", to_ugkg = FALSE)
  expect_equal(out$right_cen, c(TRUE, rep(FALSE, 9)))
  expect_equal(out$summary_recovery_status, c("recovered_summary_upper_limit", "no_upper_limit_evidence",
    "blocked_method_mismatch", rep("blocked_malformed_summary", 4), "no_summary", "no_summary", "no_upper_limit_evidence"))
  expect_equal(out$detected_val[-1], rep(10.11, 9))
})

test_that("blocked recovery preserves original censoring and invalid-code reasons", {
  x <- dplyr::bind_rows(upper_recovery_input(qualifier = "L"),
    upper_recovery_input(qualifier = "N"), upper_recovery_input(qualifier = "E"),
    upper_recovery_input(method = "AES_HF", qualifier = "G"),
    upper_recovery_input(qualifier = "G"), upper_recovery_input(value = "-10.11"),
    upper_recovery_input(value = "bad"), upper_recovery_input(value = "0"))
  out <- build_element(x, "S", to_ugkg = FALSE)
  expect_equal(out$left_cen, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(out$left_limit[c(1, 2, 6)], rep(10.11, 3))
  expect_equal(out$right_limit[c(4, 5)], c(10.11, 10))
  expect_equal(out$validation_status, c("usable", "usable", "unknown_qualifier", "usable", "usable",
    "usable", "invalid_value", "nonpositive_detection"))
  expect_equal(out$summary_recovery_status, c("blocked_censoring_conflict", "blocked_censoring_conflict",
    "blocked_reported_coding", "blocked_method_mismatch", "recovered_summary_upper_limit",
    "blocked_censoring_conflict", "blocked_reported_coding", "blocked_reported_coding"))
  expect_identical(out$qvalue_raw, x$qvalue)
})

test_that("explicit raw upper limits and disagreeing limits retain original coding", {
  x <- dplyr::bind_rows(upper_recovery_input(value = "10.1111"),
    upper_recovery_input(value = "10", qualifier = "G"),
    upper_recovery_input(value = "12", qualifier = ">"),
    upper_recovery_input(value = "9"), upper_recovery_input(summary = "0.11111, CB_IRC", value = "0.11"))
  out <- build_element(x, "S", to_ugkg = FALSE)
  expect_equal(out$right_limit[1:3], c(10, 10, 12))
  expect_true(all(out$censoring_source == "qvalue"))
  expect_equal(out$summary_recovery_status, c("original_upper_limit_retained", "original_upper_limit_retained",
    "blocked_limit_disagreement", "blocked_censoring_conflict", "blocked_summary_coding"))
  expect_equal(out$detected_val[4:5], c(9, .11))
})

test_that("recovery is unit-safe and works with an existing observational audit", {
  x <- upper_recovery_input()
  native <- build_element(x, "S", to_ugkg = FALSE)
  converted <- build_element(x, "S")
  expect_equal(converted$right_limit, native$right_limit * 1e7)
  expect_identical(converted$summary_recovery_status, native$summary_recovery_status)
  audited <- audit_summary_precision(x)
  out <- build_element(audited, "S")
  expect_identical(out[names(audited)], audited)
  expect_equal(out$right_limit, 1e8)
  empty <- build_element(x[FALSE, ], "S")
  expect_equal(nrow(empty), 0)
  expect_type(empty$summary_recovery_status, "character")
})

test_that("ordinary summaries avoid running a full source audit during parsing", {
  testthat::local_mocked_bindings(audit_summary_precision = function(...) stop("unexpected full audit"))
  x <- upper_recovery_input(summary = "10.11, CB_IRC")
  expect_equal(build_element(x, "S")$summary_recovery_status, "no_upper_limit_evidence")
  x$summary <- NULL
  expect_equal(build_element(x, "S")$summary_recovery_status, "no_summary")
})

test_that("wide comparison summaries use the same guarded recovery", {
  x <- dplyr::bind_rows(upper_recovery_input(), upper_recovery_input("AES_HF"),
    upper_recovery_input(NA_character_))
  x$lab_id <- c("001", "002", "003"); x$unit <- "ppm"
  wide <- build_metals_wide(x, "S")
  expect_equal(wide$S_mgkg, c(10, 10.11, 10.11))
})
