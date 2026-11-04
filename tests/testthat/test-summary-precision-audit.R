test_that("summary precision audit preserves data and separates source evidence", {
  x <- data.frame(id = letters[1:6], qvalue = c("0.05", "10.11", "-2", "2", "10.11", "3"),
    amethod = c("AES_HF", "CB_IRC", "AES_HF", "AES_HF", "AES_HF", "AES_HF"),
    qualifier = c("", "", "", "L", "", "<"),
    summary = c("0.054, AES_HF", "10.11111, CB_IRC", "-2, AES_HF", "-2, AES_HF",
                "10.11111, ES_SQ", "3, AES_HF"), stringsAsFactors = FALSE)
  out <- audit_summary_precision(x)
  expect_identical(out[names(x)], tibble::as_tibble(x))
  expect_equal(out$summary_precision_status, c("different_reported_value", "summary_upper_limit_evidence",
    "same_reported_value", "same_reported_value", "method_mismatch", "censoring_encoding_difference"))
  expect_equal(out$summary_numeric_delta[1:4], c(.004, -.11, 0, 0), tolerance = 1e-12)
  expect_equal(out$summary_qvalue_candidate[1:4], c(.054, 10.11111, -2, -2))
  expect_true(all(is.na(out$summary_qvalue_candidate[5:6])))
  expect_identical(x$qvalue, c("0.05", "10.11", "-2", "2", "10.11", "3"))
})

test_that("missing, malformed and invalid source information remains explicit", {
  x <- data.frame(qvalue = c("1", "1", "1", "0", "1", "1", "1"),
    amethod = c(NA, "AES_HF", "AES_HF", "AES_HF", "AES_HF", "AES_HF", "AES_HF"),
    qualifier = c("", "", "", "", "E", "", ""),
    summary = c("1, AES_HF", NA, "bad, AES_HF", ".002, AES_HF", "1, AES_HF", "1,AES_HF,junk", "1"))
  out <- audit_summary_precision(x)
  expect_equal(out$summary_precision_status, c("method_information_missing", "not_available_summary",
    "summary_coding_needs_review", "reported_coding_needs_review", "reported_coding_needs_review",
    "unparsed_summary", "unparsed_summary"))
  expect_true(all(is.na(out$summary_qvalue_candidate)))
  expect_equal(out$summary_reported_coding_status[4:5], c("nonpositive_detection", "unknown_qualifier"))
  expect_equal(audit_summary_precision(data.frame(qvalue = "1"))$summary_precision_status,
    "not_available_summary")
  trailing <- data.frame(qvalue = "1", amethod = "AES_HF", summary = "1,AES_HF,")
  expect_equal(audit_summary_precision(trailing)$summary_precision_status, "unparsed_summary")
  expect_true(is.na(audit_summary_precision(trailing)$summary_qvalue_candidate))
})

test_that("only the first same-method source entry is considered", {
  x <- data.frame(qvalue = "0.05", amethod = " aes_hf ",
    summary = "0.054, AES_HF; 20, ES_SQ")
  out <- audit_summary_precision(x)
  expect_true(out$summary_method_matches)
  expect_equal(out$summary_qvalue_candidate, .054)
  expect_equal(out$summary_first_method_raw, "AES_HF")
  wrong <- x; wrong$summary <- "20, ES_SQ; 0.054, AES_HF"
  expect_equal(audit_summary_precision(wrong)$summary_precision_status, "method_mismatch")
  expect_true(is.na(audit_summary_precision(wrong)$summary_qvalue_candidate))
})

test_that("custom columns, empty inputs and invalid schemas are handled safely", {
  x <- data.frame(v = "1", m = "AES_HF", s = "1,AES_HF", q = "")
  expect_equal(audit_summary_precision(x, "v", "m", "s", "q")$summary_precision_status,
    "same_reported_value")
  expect_equal(nrow(audit_summary_precision(x[FALSE, ], "v", "m", "s", "q")), 0L)
  expect_error(audit_summary_precision(x), "Missing reported value")
  expect_error(audit_summary_precision(x, "v", "v"), "distinct")
  expect_error(audit_summary_precision(data.frame(qvalue = I(list(1)))), "atomic")
  y <- data.frame(qvalue = "1", summary_precision_status = "user data")
  expect_error(audit_summary_precision(y), "already exist")
})

test_that("unspecified method labels cannot establish a source association", {
  x <- data.frame(qvalue = rep("0.05", 4),
    amethod = c("UNKNOWN", " unknown ", "AES_HF", NA_character_),
    summary = c("0.054, UNKNOWN", "0.054, AES_HF", "0.054, unknown", "0.054, AES_HF"))
  out <- audit_summary_precision(x)
  expect_identical(out[names(x)], tibble::as_tibble(x))
  expect_true(all(is.na(out$summary_method_matches)))
  expect_true(all(is.na(out$summary_qvalue_candidate)))
  expect_true(all(is.na(out$summary_numeric_delta)))
  expect_true(all(out$summary_precision_status == "method_information_missing"))
})
