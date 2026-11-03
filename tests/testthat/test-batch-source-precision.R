precision_batch_input <- function() {
  ids <- sprintf("%03d", 1:16)
  chemistry <- tibble::tibble(lab_id = rep(ids, 3), species = rep(c("Ag", "Mn", "S"), each = 16),
    unit = rep(c("ppm", "pct", "pct"), each = 16),
    qvalue = c(as.character(1:16), as.character(seq(.05, .20, by = .01)), as.character(1:16)),
    amethod = rep(c("AES_HF", "AES_HF", "CB_IRC"), each = 16), qualifier = "")
  chemistry$summary <- paste(chemistry$qvalue, chemistry$amethod, sep = ", ")
  chemistry$lab_id[chemistry$lab_id == "002"] <- " 002 "
  chemistry$qvalue[3] <- "-3"; chemistry$summary[3] <- "-3, AES_HF"
  chemistry$summary[2] <- "2.4, AES_HF"
  chemistry$summary[17:21] <- paste(c("0.054", "0.064", "0.074", "0.084", "0.094"), "AES_HF", sep = ", ")
  chemistry$amethod[22] <- "UNKNOWN"; chemistry$summary[22] <- "0.104, UNKNOWN"
  chemistry$qvalue[33] <- "10.11"; chemistry$summary[33] <- "10.11111, CB_IRC"
  chemistry$summary[48] <- NA_character_
  chemistry$unit[16] <- "unsupported"
  chemistry <- chemistry[as.vector(rbind(33:48, 1:16, 17:32)), ]
  metadata <- tibble::tibble(lab_id = ids, state = c(rep("WI", 8), rep("IL", 7), NA_character_),
    x_m = c(NA_real_, seq(1000, 15000, by = 1000)), y_m = 0)
  list(chemistry = chemistry, metadata = metadata)
}

test_that("source review is opt-in and preserves every existing analysis component", {
  input <- precision_batch_input(); saved <- input
  for (stratify in c(FALSE, TRUE)) for (overlap in c(FALSE, TRUE)) {
    args <- list(chem_df = input$chemistry, main_df = input$metadata, group_cols = "state",
      cell_sizes = c(1000, 5000), x_col = "x_m", y_col = "y_m", n_offsets = 1,
      stratify_methods = stratify, audit_method_overlap = overlap)
    baseline <- do.call(analyze_elements, args)
    result <- do.call(analyze_elements, c(args, list(audit_source_precision = TRUE)))
    expect_named(result, c(names(baseline), "source_precision", "source_precision_summary"))
    expect_equal(result[names(baseline)], baseline)
    expect_true(all(is.na(result$declustering$ci_low)))
  }
  expect_identical(input, saved)
})

test_that("row audit retains original order, raw IDs, native units and excluded inputs", {
  input <- precision_batch_input()
  result <- analyze_elements(input$chemistry, main_df = input$metadata,
    group_cols = "state", audit_source_precision = TRUE)
  rows <- result$source_precision
  expect_identical(rows[names(input$chemistry)], input$chemistry)
  standalone <- audit_summary_precision(input$chemistry)
  expect_identical(rows[names(standalone)], standalone)
  expect_identical(rows$state, input$metadata$state[match(trimws(rows$lab_id), input$metadata$lab_id)])
  expect_false(any(c("x_m", "y_m") %in% names(rows)))
  expect_true(any(rows$lab_id == " 002 "))
  mn <- rows[rows$species == "Mn" & trimws(rows$lab_id) == "001", ]
  ag <- rows[rows$species == "Ag" & trimws(rows$lab_id) == "002", ]
  expect_equal(mn$summary_numeric_delta, .004, tolerance = 1e-12)
  expect_equal(ag$summary_numeric_delta, .4, tolerance = 1e-12)
  expect_equal(mn$unit, "pct")
  expect_equal(ag$unit, "ppm")
  expect_equal(sum(rows$unit == "unsupported"), 1)
  expect_equal(sum(result$estimates$n_excluded), 1)
  expect_equal(sum(result$source_precision_summary$n_observations), 48)
  expect_equal(table(rows$summary_precision_status), table(c(rep("same_reported_value", 39),
    rep("different_reported_value", 6), "method_information_missing",
    "summary_upper_limit_evidence", "not_available_summary")))
})

test_that("review compares raw upper-limit coding without changing parser recovery", {
  input <- precision_batch_input()
  result <- analyze_elements(input$chemistry, audit_source_precision = TRUE)
  row <- result$source_precision[result$source_precision$species == "S" &
    trimws(result$source_precision$lab_id) == "001", ]
  expect_equal(row$qvalue, "10.11")
  expect_equal(row$summary_precision_status, "summary_upper_limit_evidence")
  expect_equal(row$summary_numeric_delta, -.11, tolerance = 1e-12)
  expect_equal(row$summary_qvalue_candidate, 10.11111)
  parsed <- build_element(input$chemistry, "S")
  parsed <- parsed[trimws(parsed$lab_id) == "001", ]
  expect_equal(parsed$censoring_source, "summary_upper_limit")
  expect_equal(parsed$right_limit, 10 * 1e7)
  expect_equal(parsed$qvalue_raw, "10.11")
})

test_that("status counts follow requested elements and metadata domains including NA", {
  input <- precision_batch_input()
  result <- analyze_elements(input$chemistry, elements = c("Mn", "Au"),
    main_df = input$metadata, group_cols = "state", audit_source_precision = TRUE)
  rows <- result$source_precision; counts <- result$source_precision_summary
  expect_equal(nrow(rows), 16)
  expect_true(all(rows$species == "Mn"))
  expect_false(any(counts$element == "Au"))
  expect_true(anyNA(counts$state))
  for (i in seq_len(nrow(counts))) {
    domain <- if (is.na(counts$state[i])) is.na(rows$state) else !is.na(rows$state) & rows$state == counts$state[i]
    expect_equal(counts$n_observations[i], sum(domain & rows$summary_precision_status == counts$summary_precision_status[i]))
  }
  absent <- analyze_elements(input$chemistry, elements = "Au", audit_source_precision = TRUE)
  expect_equal(nrow(absent$source_precision), 0)
  expect_equal(nrow(absent$source_precision_summary), 0)
  expect_true(all(c("element", "summary_precision_status", "n_observations") %in% names(absent$source_precision_summary)))
  missing <- input$chemistry; missing$summary <- NULL
  unknown <- analyze_elements(missing, audit_source_precision = TRUE)
  expect_true(all(unknown$source_precision$summary_precision_status == "not_available_summary"))
  expect_true(all(is.na(unknown$source_precision$summary_qvalue_candidate)))
})

test_that("source audit rejects invalid flags and ambiguous output schemas", {
  input <- precision_batch_input()
  for (bad in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
    expect_error(analyze_elements(input$chemistry, audit_source_precision = bad),
      "audit_source_precision must be TRUE or FALSE")
  }
  for (column in c("summary_precision_status", "summary_numeric_delta", "n_observations")) {
    metadata <- input$metadata; metadata[[column]] <- "domain"
    expect_error(analyze_elements(input$chemistry, main_df = metadata, group_cols = column,
      audit_source_precision = TRUE), "Grouping columns conflict with source-precision fields")
  }
  already <- audit_summary_precision(input$chemistry)
  expect_error(analyze_elements(already, audit_source_precision = TRUE), "already exist")
  duplicate <- dplyr::bind_rows(input$chemistry, input$chemistry[1, ])
  expect_error(analyze_elements(duplicate, audit_source_precision = TRUE), "duplicate lab_id")
})
