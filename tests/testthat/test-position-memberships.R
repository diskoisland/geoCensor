position_metadata <- function() tibble::tibble(
  lab_id = c("001", " 002 ", "003", "004", "005", "006", "007"),
  state = c("A", "B", "A", NA_character_, "A", "B", "B"),
  depth = c(" 0-10 ", "10-20", NA_character_, "0-10", "0-10", "0-10", "0-10"),
  x = c(1, 1, 2, NA, Inf, 3, 3), y = c(1, 1, 2, 2, 1, 3, 3))

position_audit <- function(x = position_metadata(), ...) sampling_support_audit(x,
  group_cols = "state", support_cols = "depth", x_col = "x", y_col = "y", ...)

test_that("membership opt-in preserves original metadata and the four default components", {
  x <- position_metadata(); saved <- x
  base <- position_audit(x)
  out <- position_audit(x, include_memberships = TRUE)
  expect_named(out, c(names(base), "memberships", "coordinate_groups"))
  expect_identical(out[names(base)], base)
  expect_identical(out$memberships[names(x)], x)
  expect_identical(x, saved)
  expect_equal(out$memberships$n_records_at_position, c(2L, 2L, 1L, NA, NA, 2L, 2L))
  expect_equal(out$memberships$position_membership_status, c(rep("recorded_position_group", 3),
    rep("unavailable_coordinates", 2), rep("recorded_position_group", 2)))
  expect_equal(out$memberships$review_shared_position, c(TRUE, TRUE, NA, NA, NA, TRUE, TRUE))
  expect_true(all(is.na(out$memberships$recorded_position_group[4:5])))
  expect_equal(nrow(out$coordinate_groups), 3)
})

test_that("global groups retain cross-domain positions without changing domain summaries", {
  out <- position_audit(include_memberships = TRUE)
  groups <- out$coordinate_groups
  shared <- groups[groups$x == 1, ]
  expect_equal(shared$n_records_at_position, 2)
  expect_equal(shared$n_analysis_domains, 2)
  expect_true(shared$mixed_recorded_support)
  expect_false(shared$has_support_gaps)
  singleton <- groups[groups$x == 2, ]
  expect_equal(singleton$n_records_at_position, 1)
  expect_true(singleton$has_support_gaps)
  expect_true(is.na(singleton$review_shared_position))
  expect_equal(nrow(out$locations[out$locations$x == 1, ]), 2)
  expect_true(all(out$locations$n_samples[out$locations$x == 1] == 1))
  expect_identical(out$memberships$recorded_position_group[1], out$memberships$recorded_position_group[2])
  ungrouped <- sampling_support_audit(position_metadata(), support_cols = "depth", x_col = "x", y_col = "y",
    include_memberships = TRUE)
  expect_true(all(ungrouped$coordinate_groups$n_analysis_domains == 1L))
  x <- position_metadata(); x$state[2] <- NA_character_
  expect_equal(position_audit(x, include_memberships = TRUE)$coordinate_groups$n_analysis_domains[1], 2)
})

test_that("exact numeric grouping is stable under reordering and preserves close coordinates", {
  x <- position_metadata()
  first <- position_audit(x, include_memberships = TRUE)
  reordered <- x[c(7, 4, 1, 6, 3, 5, 2), ]
  second <- position_audit(reordered, include_memberships = TRUE)
  expect_identical(second$coordinate_groups, first$coordinate_groups)
  expect_identical(second$memberships$lab_id, reordered$lab_id)
  expect_identical(second$memberships$recorded_position_group,
    first$memberships$recorded_position_group[match(reordered$lab_id, x$lab_id)])
  for (values in list(c(2^53 - 1, 2^53, 2^53 - 1), c(1, 1 + .Machine$double.eps, 1))) {
    z <- tibble::tibble(lab_id = c("001", "002", "003"), depth = "A", x = values, y = 0)
    result <- sampling_support_audit(z, x_col = "x", y_col = "y", include_memberships = TRUE)
    expect_equal(nrow(result$coordinate_groups), 2)
    expect_identical(result$memberships$recorded_position_group[1], result$memberships$recorded_position_group[3])
    expect_false(result$memberships$recorded_position_group[1] == result$memberships$recorded_position_group[2])
  }
})

test_that("empty and unavailable coordinates never become an invented shared group", {
  x <- position_metadata(); x$x <- NA_real_; x$y <- Inf
  out <- position_audit(x, include_memberships = TRUE)
  expect_equal(nrow(out$memberships), 7)
  expect_equal(nrow(out$coordinate_groups), 0)
  expect_true(all(is.na(out$memberships$recorded_position_group)))
  expect_true(all(is.na(out$memberships$n_records_at_position)))
  expect_true(all(is.na(out$memberships$review_shared_position)))
  expect_true(all(out$memberships$position_membership_status == "unavailable_coordinates"))
  empty <- position_audit(x[FALSE, ], include_memberships = TRUE)
  expect_equal(nrow(empty$memberships), 0)
  expect_equal(nrow(empty$coordinate_groups), 0)
  expect_type(empty$memberships$recorded_position_group, "character")
  expect_type(empty$memberships$n_records_at_position, "integer")
})

test_that("optional output fields and temporary indexing cannot overwrite source columns", {
  x <- position_metadata()
  for (name in c("recorded_position_group", "n_records_at_position", "position_membership_status",
                 "review_shared_position", "n_analysis_domains")) {
    original <- x; original[[name]] <- "source evidence"
    expect_error(position_audit(original, include_memberships = TRUE), "already exist")
  }
  x$..geocensor_position_index <- "source grouping"
  out <- sampling_support_audit(x, group_cols = "..geocensor_position_index", support_cols = "depth",
    x_col = "x", y_col = "y", include_memberships = TRUE)
  expect_identical(out$memberships[names(x)], x)
  expect_true(all(out$coordinate_groups$n_analysis_domains == 1L))
  names(x)[1] <- "sample_key"
  custom <- sampling_support_audit(x, id_col = "sample_key", x_col = "x", y_col = "y", include_memberships = TRUE)
  expect_identical(custom$memberships$sample_key, x$sample_key)
})

test_that("membership flags and coordinate types are validated explicitly", {
  for (bad in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
    expect_error(position_audit(include_memberships = bad), "include_memberships must be TRUE or FALSE")
  }
  expect_error(sampling_support_audit(position_metadata(), include_memberships = TRUE), "requires both coordinate")
  for (values in list(as.complex(1:7), matrix(1:7, ncol = 1), as.character(1:7))) {
    x <- position_metadata(); x$x <- values
    expect_error(position_audit(x, include_memberships = TRUE), "distinct numeric")
  }
})
