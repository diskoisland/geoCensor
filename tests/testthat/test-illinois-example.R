test_that("Illinois example selects by metadata and preserves censored and source records", {
  root <- tempfile("illinois_example_"); dir.create(root)
  main <- file.path(root, "main.csv"); chemistry <- file.path(root, "bestvalue.csv")
  projected <- file.path(root, "projected.csv"); dest <- file.path(root, "result")
  ids <- sprintf("%03d", seq_len(12))
  metadata <- data.frame(lab_id = c(ids, "999"), state = "IL",
    project_name = c(rep("NAWQA Upper Illinois R. Basin", 12), "Other project"),
    long_wgs84 = c(rep("-89", 6), rep("-88", 6), "-90"), lat_wgs84 = "41",
    addl_attr = "unverified sediment fraction", date_collect = NA_character_)
  readr::write_csv(metadata, main)
  readr::write_csv(data.frame(lab_id = ids, x_m = seq_len(12)*1000,
    y_m = rep(c(0, 2000), 6)), projected)
  chem <- data.frame(lab_id = c(ids, ids, "999"), species = c(rep("Ag", 12), rep("Au", 12), "Ag"),
    unit = "ppm", qvalue = c("-2", as.character(2:12), rep("-5", 12), "9999"),
    amethod = "MS_HF", qualifier = "", source_record_number = "000123")
  chem$qvalue[12] <- "4.11110"
  readr::write_csv(chem, chemistry)
  run <- run_analysis_example("illinois_case_study.R", c(main, chemistry, dest,
    paste0("--projected-metadata=", projected), "--chunk-size=5"))
  expect_equal(run$status, 0L, info = paste(run$log, collapse = "\n"))
  if (run$status != 0L) return(invisible(NULL))
  result <- readRDS(file.path(dest, "analysis.rds"))
  input <- read_usgs_chemistry(file.path(dest, "selected_chemistry.csv"))
  expect_identical(input$lab_id, chem$lab_id[1:24])
  expect_identical(input$qvalue[12], "4.11110")
  expect_true(all(input$source_record_number == "000123"))
  expect_setequal(result$estimates$element, c("Ag", "Au"))
  expect_true(is.na(result$estimates$geomean[result$estimates$element == "Au"]))
  expect_true(all(is.na(result$declustering$ci_low)))
  expect_true(all(is.na(result$declustering$ci_high)))
  bounds <- readRDS(file.path(dest, "Ag_cdf_bounds.rds"))
  expect_equal(bounds$support$n_input, 12)
  expect_true(all(bounds$cdf$cdf_lower <= bounds$cdf$cdf_upper))
  expect_true(any(bounds$cdf$censoring_gap > 0))
  inventory <- readr::read_csv(file.path(dest, "inventory.csv"), show_col_types = FALSE)
  expect_equal(inventory$value[inventory$item == "metadata_records"], 12)
  expect_equal(inventory$value[inventory$item == "chemistry_records"], 24)
  expect_equal(inventory$value[inventory$item == "recorded_geographic_positions"], 2)
  before <- tools::md5sum(list.files(dest, full.names = TRUE))
  rerun <- run_analysis_example("illinois_case_study.R", c(main, chemistry, dest,
    paste0("--projected-metadata=", projected)))
  expect_true(rerun$status != 0L)
  expect_identical(tools::md5sum(list.files(dest, full.names = TRUE)), before)
  # Malformed unselected records must still abort before output creation.
  cat("999,Pb,ppm,1,MS_HF,,000123,extra\n", file = chemistry, append = TRUE)
  bad <- file.path(root, "bad")
  run <- run_analysis_example("illinois_case_study.R", c(main, chemistry, bad,
    paste0("--projected-metadata=", projected), "--chunk-size=5"))
  expect_true(run$status != 0L)
  expect_false(file.exists(bad))
})
