multi_chemistry <- function() {
  tibble::tibble(
    lab_id = rep(sprintf("%03d", 1:8), 3),
    species = rep(c("Ag", "Pb", "Zn"), each = 8),
    qvalue = rep(c("-0.02", "0.03", "0.05", "0.08", "0.1", "0.15", "0.2", "0.4"), 3),
    unit = rep(c("ppm", "ppb", "unknown"), each = 8)
  )
}

test_that("all species are analyzed with explicit unsupported-unit results", {
  x <- multi_chemistry()
  out <- analyze_elements(x)
  expect_setequal(out$estimates$element, c("Ag", "Pb", "Zn"))
  expect_equal(nrow(out$distributions), 9)
  for (element in c("Ag", "Pb")) {
    reference <- fit_cen_lognorm(build_element(x, element))
    expect_equal(out$estimates$geomean[out$estimates$element == element], reference$geomean)
  }
  expect_equal(out$estimates$n_excluded[out$estimates$element == "Zn"], 8)
  expect_equal(out$estimates$method[out$estimates$element == "Zn"], "insufficient_data")
  expect_equal(out$validation$n_observations[out$validation$validation_status == "unsupported_unit"], 8)
  expect_equal(nrow(out$declustering), 0)
})

test_that("requested absent elements receive an insufficient-data row", {
  out <- analyze_elements(multi_chemistry(), elements = c("Ag", "Au"))
  expect_equal(out$estimates$n_input, c(8, 0))
  expect_equal(out$estimates$method[2], "insufficient_data")
  expect_setequal(out$validation$element, "Ag")
})

test_that("declustering weights are recalculated per element and domain", {
  x <- multi_chemistry()
  metadata <- tibble::tibble(lab_id = sprintf("%03d", 1:8),
                            state = "WI", x_m = c(rep(0, 6), 100, 200), y_m = 0)
  out <- analyze_elements(x, main_df = metadata, group_cols = "state",
                          cell_sizes = c(20, 40), x_col = "x_m", y_col = "y_m", n_offsets = 2)
  expect_equal(nrow(out$declustering), 6)
  expect_true(all(out$estimates$state == "WI"))
  expect_true(all(out$joins$n_matched == 8))
  expect_equal(out$declustering$n_total[out$declustering$element == "Zn"], c(0, 0))
  expect_true(all(is.na(out$declustering$weighted_exact_fraction[out$declustering$element == "Zn"])))
  valid <- out$declustering$element %in% c("Ag", "Pb")
  expect_equal(out$declustering$exact_fraction[valid], rep(7 / 8, 4))
  expect_equal(out$declustering$weighted_exact_fraction[valid], rep(17 / 18, 4))
  expect_equal(out$declustering$weighted_left_fraction[valid], rep(1 / 18, 4))
  expect_equal(out$declustering$weighted_right_fraction[valid], rep(0, 4))
  for (element in c("Ag", "Pb")) {
    data <- join_soil_main(build_element(x, element), metadata, extra_cols = c("x_m", "y_m"))
    ref <- declustering_sensitivity(data, c(20, 40), "x_m", "y_m", n_offsets = 2)
    expect_equal(out$declustering$geomean[out$declustering$element == element], ref$geomean)
  }
})

test_that("batch analyses enforce sample identity and grouping requirements", {
  x <- multi_chemistry()
  x$lab_id[2] <- x$lab_id[1]
  expect_error(analyze_elements(x), "Ag: duplicate")
  x$lab_id[2] <- NA_character_
  expect_error(analyze_elements(x), "missing or blank")
  expect_error(analyze_elements(multi_chemistry(), group_cols = "medium"), "Missing grouping")
  expect_error(analyze_elements(multi_chemistry(), cell_sizes = 100), "projected")
})

test_that("separate domains retain their own baseline and sample counts", {
  first <- multi_chemistry()[1:16, ]
  second <- first
  second$lab_id <- paste0(second$lab_id, "b")
  second$qvalue <- as.character(seq_len(nrow(second)))
  x <- dplyr::bind_rows(first, second)
  metadata <- tibble::tibble(
    lab_id = c(sprintf("%03d", 1:8), paste0(sprintf("%03d", 1:8), "b")),
    domain = rep(c("first", "second"), each = 8),
    x_m = c(rep(0, 6), 100, 200, seq(100, 800, by = 100)), y_m = 0
  )
  out <- analyze_elements(x, main_df = metadata, group_cols = "domain",
                          cell_sizes = 20, x_col = "x_m", y_col = "y_m", n_offsets = 2)
  expect_equal(nrow(out$estimates), 4)
  expect_true(all(out$estimates$n_total == 8))
  for (element in c("Ag", "Pb")) {
    for (domain in c("first", "second")) {
      ids <- metadata$lab_id[metadata$domain == domain]
      observations <- join_soil_main(build_element(x[x$lab_id %in% ids, ], element),
                                      metadata, extra_cols = c("x_m", "y_m"))
      reference <- declustering_sensitivity(observations, 20, "x_m", "y_m", n_offsets = 2)
      matched <- out$declustering$element == element & out$declustering$domain == domain
      expect_equal(out$declustering$geomean[matched], reference$geomean)
    }
  }
})
