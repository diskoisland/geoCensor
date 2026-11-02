test_that("build_ag_from_chem handles missing qualifier and converts units", {
  chem <- tibble::tibble(
    lab_id = c("a", "b", "c", "d"),
    species = c("Ag", "Ag", "Ag", "Pb"),
    qvalue = c("-2", "4", "0.01", "10"),
    unit = c("ppm", "ppm", "pct", "ppm")
  )

  ag <- build_ag_from_chem(chem, element = "Ag", to_ugkg = TRUE)

  expect_equal(nrow(ag), 3)
  expect_true("qualifier" %in% names(ag))
  expect_true(ag$left_cen[ag$lab_id == "a"])
  expect_equal(ag$left_limit[ag$lab_id == "a"], 2000)
  expect_equal(ag$detected_val[ag$lab_id == "b"], 4000)
  expect_equal(ag$detected_val[ag$lab_id == "c"], 100000)
})

test_that("build_metals_wide substitutes censored values and handles all-NA groups", {
  chem <- tibble::tibble(
    lab_id = c("a", "a", "b", "c"),
    species = c("Ag", "Pb", "Ag", "Ag"),
    qvalue = c("-2", "10", "0.01", "5"),
    unit = c("ppm", "ppm", "pct", "ppb")
  )

  wide <- build_metals_wide(chem, elements = c("Ag", "Pb"))

  expect_equal(wide$Ag_mgkg[wide$lab_id == "a"], 1)
  expect_equal(wide$Pb_mgkg[wide$lab_id == "a"], 10)
  expect_equal(wide$Ag_mgkg[wide$lab_id == "b"], 100)
  expect_equal(wide$Ag_mgkg[wide$lab_id == "c"], 0.005)
})

test_that("build_ag_metals_state joins metadata and filters state names", {
  chem <- tibble::tibble(
    lab_id = c("a", "b"),
    species = c("Ag", "Ag"),
    qvalue = c("1", "2"),
    unit = c("ppm", "ppm")
  )
  main <- tibble::tibble(
    lab_id = c("a", "b"),
    state = c("WI", "ID"),
    lat_wgs84 = c(44, 45),
    long_wgs84 = c(-89, -116)
  )

  out <- build_ag_metals_state(chem, main, elements = "Ag", state = "Wisconsin")

  expect_equal(out$lab_id, "a")
  expect_equal(out$state, "WI")
  expect_equal(out$Ag_mgkg, 1)
})
