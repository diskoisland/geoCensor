test_that("norm_state normalizes names, abbreviations, and blanks", {
  expect_equal(
    norm_state(c("Wisconsin", "wi", "", "Not a state")),
    c("WI", "WI", NA_character_, NA_character_)
  )
})

test_that("read_clean_main normalizes state and invalid coordinates", {
  main <- tibble::tibble(
    lab_id = c("a", "b"),
    state = c("Wisconsin", "XX"),
    country = c(" US ", "USA"),
    lat_wgs84 = c(44, 100),
    long_wgs84 = c(-89, -200),
    date_collect = c("01/02/2020", "bad"),
    date_submitted = c("02/03/2020", "bad"),
    depth = c("0-10 cm", "5 m")
  )
  path <- tempfile(fileext = ".csv")
  readr::write_csv(main, path)

  out <- read_clean_main(path)

  expect_equal(out$state, c("WI", NA_character_))
  expect_equal(out$country, c("US", "USA"))
  expect_equal(out$lat_wgs84, c(44, NA_real_))
  expect_equal(out$long_wgs84, c(-89, NA_real_))
  expect_equal(out$depth_num, c(10, 5))
  expect_equal(out$depth_unit, c("cm", "m"))
})

test_that("join_soil_main and join_sed_main retain requested metadata", {
  elem <- tibble::tibble(lab_id = c("a", "b"), detected_val = c(1, 2))
  main <- tibble::tibble(
    lab_id = c("a", "b"),
    state = c("WI", "ID"),
    country = c("US", "US"),
    extra = c("x", "y")
  )

  soil <- join_soil_main(elem, main, extra_cols = c("state", "extra"))
  sed <- join_sed_main(elem, main, extra_cols = c("country"))

  expect_equal(names(soil), c("lab_id", "detected_val", "state", "extra"))
  expect_equal(names(sed), c("lab_id", "detected_val", "country"))
})
