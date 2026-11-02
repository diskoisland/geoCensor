test_that("spatial_coverage_by_state counts valid coordinates and grid cells", {
  df <- tibble::tibble(
    state = c("WI", "WI", "WI", "ID", "ID", "ID", "ID", NA, "TX", "TX"),
    lat_wgs84 = c(43.01, 43.02, 43.76, 44.1, 44.9, 45.2, 45.3, 40, 95, 31),
    long_wgs84 = c(-89.41, -89.42, -90.1, -113.1, -113.2, -114.0, -114.1, -90, -100, -190)
  )

  out <- spatial_coverage_by_state(df, grid_deg = 0.5)

  expect_equal(out$state, c("ID", "WI"))
  expect_equal(out$n_spatial[out$state == "WI"], 3)
  expect_equal(out$occupied_grid_cells[out$state == "WI"], 2)
  expect_equal(out$n_spatial[out$state == "ID"], 4)
  expect_equal(out$occupied_grid_cells[out$state == "ID"], 4)
  expect_equal(out$grid_deg, c(0.5, 0.5))
  expect_true(out$bbox_area_deg2[out$state == "ID"] > out$bbox_area_deg2[out$state == "WI"])
})

test_that("spatial_coverage_by_state supports custom column names", {
  df <- tibble::tibble(
    st = c("CA", "CA"),
    y = c(36.1, 37.1),
    x = c(-119.2, -119.8)
  )

  out <- spatial_coverage_by_state(
    df,
    state_col = "st",
    lat_col = "y",
    lon_col = "x",
    grid_deg = 1
  )

  expect_equal(out$state, "CA")
  expect_equal(out$n_spatial, 2)
  expect_equal(out$occupied_grid_cells, 2)
})

test_that("spatial_coverage_by_state handles empty valid coordinates", {
  df <- tibble::tibble(
    state = c("WI", NA),
    lat_wgs84 = c(NA_real_, 95),
    long_wgs84 = c(-89, -200)
  )

  out <- spatial_coverage_by_state(df)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0)
  expect_named(out, c(
    "state", "n_spatial", "grid_deg", "occupied_grid_cells",
    "samples_per_grid_cell", "lat_min", "lat_max", "lon_min", "lon_max",
    "lat_range", "lon_range", "bbox_area_deg2"
  ))
})

test_that("spatial_coverage_by_state validates inputs", {
  df <- tibble::tibble(state = "WI", lat_wgs84 = 43)

  expect_error(
    spatial_coverage_by_state(df),
    "Missing columns: long_wgs84"
  )
  expect_error(
    spatial_coverage_by_state(
      tibble::tibble(state = "WI", lat_wgs84 = 43, long_wgs84 = -89),
      grid_deg = 0
    ),
    "grid_deg must be one positive number"
  )
})

test_that("spatial coverage reads coordinate labels rather than factor level codes", {
  x <- tibble::tibble(state = c(" WI ", "WI", " "),
                      lat_wgs84 = factor(c("43.1", "43.9", "40")),
                      long_wgs84 = factor(c("-89.1", "-89.9", "-100")))
  out <- spatial_coverage_by_state(x)
  expect_equal(out$state, "WI")
  expect_equal(out$n_spatial, 2)
  expect_equal(out$lat_min, 43.1)
  expect_equal(out$lat_max, 43.9)
  expect_equal(out$lon_min, -89.9)
  expect_equal(out$lon_max, -89.1)
})
