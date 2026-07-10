test_that("plot_ranked_bar returns a ggplot object", {
  df <- tibble::tibble(
    state = c("WI", "ID"),
    geomean = c(10, 20),
    ci_low = c(8, 15),
    ci_high = c(12, 25)
  )

  expect_s3_class(plot_ranked_bar(df), "ggplot")
})

test_that("plot_us_map supports custom element and medium labels", {
  df <- tibble::tibble(
    state = c("WI", "ID"),
    geomean = c(10, 20),
    ci_low = c(8, 15),
    ci_high = c(12, 25)
  )

  p <- plot_us_map(df, element_label = "Ag", medium_label = "sediments")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "State-level Ag geometric mean in U.S. sediments (MLE)")
  expect_equal(p$scales$scales[[1]]$name, "Ag (ug kg^-1)")
})

test_that("state_hotspot_map returns a ggplot object for matching points", {
  skip_if_not_installed("maps")

  df <- tibble::tibble(
    state = c("WI", "WI"),
    long_wgs84 = c(-89, -90),
    lat_wgs84 = c(44, 45),
    Ag_mgkg = c(2.5, 3)
  )

  expect_s3_class(state_hotspot_map(df, "Wisconsin", thr_min = 2), "ggplot")
})
