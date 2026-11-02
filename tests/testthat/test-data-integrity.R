test_that("invalid encodings carry reasons and cannot enter estimates", {
  x <- parse_usgs_q(
    c("-2", "4.1111", "0", "bad", NA, "0", "5", "2", "Inf"),
    c("G", "L", "", "", "", "L", "?", "", "")
  )
  expect_equal(x$validation_status, c(
    "conflicting_qualifiers", "conflicting_qualifiers", "nonpositive_detection",
    "invalid_value", "missing_value", "invalid_limit", "unknown_qualifier",
    "usable", "invalid_value"
  ))
  expect_true(all(is.na(x$detected_val[-8])))
  expect_equal(fit_cen_lognorm(x)$n_total, 1)
})

test_that("equivalent units convert detections and both censoring limits equally", {
  units <- c("ppm", " mg/kg ", "ppb", "ug/kg", "\u00b5g/kg", "\u03bcg/kg", "pct", "%", "percent")
  factors <- c(1000, 1000, 1, 1, 1, 1, 1e7, 1e7, 1e7)
  for (value in c("2", "-2", "2.1111")) {
    chem <- tibble::tibble(lab_id = seq_along(units), species = "Ag", unit = units, qvalue = value)
    out <- build_element(chem, "Ag")
    column <- switch(value, "2" = "detected_val", "-2" = "left_limit", "2.1111" = "right_limit")
    expect_equal(out[[column]], 2 * factors)
    expect_true(all(out$validation_status == "usable"))
    wide <- build_metals_wide(chem, "Ag", lod_factor = 1)
    expect_equal(wide$Ag_mgkg, out[[column]] / 1000)
  }
})

test_that("unsupported units are explicit and original units can be retained", {
  chem <- tibble::tibble(lab_id = "001", species = "Ag", unit = "mg/L", qvalue = "2")
  out <- build_element(chem, "Ag")
  expect_equal(out$validation_status, "unsupported_unit")
  expect_true(is.na(out$detected_val))
  expect_equal(build_element(chem, "Ag", to_ugkg = FALSE)$detected_val, 2)
  expect_error(build_element(data.frame(), "Ag"), "Missing chemistry columns")
  expect_error(build_element(chem, c("Ag", "Pb")), "one nonempty")
  expect_warning(wide <- build_metals_wide(chem, "Ag"), "chemistry rows excluded")
  expect_equal(attr(wide, "validation_summary")$validation_status, "unsupported_unit")
  expect_true(is.na(wide$Ag_mgkg))
})

test_that("CSV readers preserve sample identifiers with leading zeroes", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("lab_id,species,qvalue,unit", "001,Ag,2,ppm", "010,Ag,3,ppm"), path)
  expect_equal(build_ag_from_csv(path)$lab_id, c("001", "010"))
})

test_that("metadata cannot multiply chemistry rows", {
  chemistry <- tibble::tibble(lab_id = c("a", "b"), detected_val = c(1, 2))
  metadata <- tibble::tibble(lab_id = c("a", "a"), state = c("WI", "ID"))
  expect_error(join_soil_main(chemistry, metadata), "Duplicate metadata")
  expect_error(join_sed_main(chemistry, metadata), "Duplicate metadata")
  metadata$lab_id[2] <- " "
  expect_error(join_soil_main(chemistry, metadata), "missing or blank")
  metadata$lab_id <- c("a", "b")
  out <- join_soil_main(chemistry, metadata)
  expect_equal(nrow(out), nrow(chemistry))
  expect_equal(attr(out, "join_audit")$n_unmatched, 0)
})

test_that("unmatched and repeated chemistry IDs are retained and reported", {
  chemistry <- tibble::tibble(lab_id = c("a", "a", "missing", NA_character_))
  metadata <- tibble::tibble(lab_id = "a", state = "WI")
  expect_warning(expect_warning(
    out <- join_soil_main(chemistry, metadata), "no matching metadata"
  ), "repeated chemistry")
  expect_equal(nrow(out), 4)
  expect_equal(out$state, c("WI", "WI", NA_character_, NA_character_))
  expect_equal(attr(out, "join_audit")$n_unmatched, 2)
  expect_equal(attr(out, "join_audit")$n_repeated_chemistry, 1)
})

test_that("wide summaries reject duplicate results instead of averaging", {
  chem <- tibble::tibble(lab_id = c("a", "a"), species = "Ag", unit = "ppm", qvalue = c(1, 100))
  expect_error(build_metals_wide(chem, "Ag"), "Duplicate lab_id/species")
  expect_error(build_metals_wide(chem[1, ], "Ag", lod_factor = NA_real_), "lod_factor")
})

test_that("rounded best-value exports recover explicit upper limits from summary", {
  chem <- tibble::tibble(
    lab_id = c("a", "b", "c", "d"), species = "Ag", unit = "ppm",
    qvalue = c("10.11", "10.11", "10.11", "20.11"),
    summary = c("10.11111, MS_HF", "10.11, MS_HF", NA, "10.11111, MS_HF")
  )
  out <- build_element(chem, "Ag")
  expect_equal(out$right_cen, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(out$right_limit[1], 10000)
  expect_equal(out$qvalue_raw[1], "10.11")
  expect_equal(out$censoring_source[1], "summary_upper_limit")
  expect_equal(out$detected_val[2:4], c(10110, 10110, 20110))
  expect_equal(build_metals_wide(chem, "Ag")$Ag_mgkg, c(10, 10.11, 10.11, 20.11))
  expect_equal(parse_usgs_q(c("10.1111", "10.11111"))$right_limit, c(10, 10))
  expect_true(parse_usgs_q(2, "N")$left_cen)
})
