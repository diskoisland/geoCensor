test_that("mle_by_state fits by state and preserves insufficient data rows", {
  enough <- tibble::tibble(
    state = "WI",
    country = "US",
    detected_val = c(NA, NA, 80, 100, 120, 150, 200, NA),
    left_cen = c(TRUE, TRUE, rep(FALSE, 6)),
    left_limit = c(40, 60, rep(NA_real_, 6)),
    right_cen = c(rep(FALSE, 7), TRUE),
    right_limit = c(rep(NA_real_, 7), 300)
  )
  sparse <- tibble::tibble(
    state = "ID",
    country = "US",
    detected_val = c(10, 20),
    left_cen = c(FALSE, FALSE),
    left_limit = c(NA_real_, NA_real_),
    right_cen = c(FALSE, FALSE),
    right_limit = c(NA_real_, NA_real_)
  )

  out <- mle_by_state(dplyr::bind_rows(enough, sparse))

  expect_setequal(out$state, c("WI", "ID"))
  expect_true(is.finite(out$geomean[out$state == "WI"]))
  expect_true(is.na(out$geomean[out$state == "ID"]))
})

test_that("add_substitution adds ug/kg and mg/kg substitution values", {
  df <- tibble::tibble(
    detected_val = c(100, NA, NA),
    left_cen = c(FALSE, TRUE, FALSE),
    left_limit = c(NA, 40, NA),
    right_cen = c(FALSE, FALSE, TRUE),
    right_limit = c(NA, NA, 300)
  )

  out <- add_substitution(df)

  expect_equal(out$val_sub_ugkg, c(100, 20, 300))
  expect_equal(out$val_sub_mgkg, c(0.1, 0.02, 0.3))
})

test_that("robust_cor_table computes bootstrap Spearman correlations", {
  set.seed(1)
  df <- tibble::tibble(
    Ag_mgkg = 1:20,
    Pb_mgkg = 1:20,
    Zn_mgkg = c(rep(NA_real_, 15), 1:5)
  )

  out <- robust_cor_table(
    df,
    target = "Ag_mgkg",
    metals = c("Pb_mgkg", "Zn_mgkg"),
    n_boot = 25,
    min_pairs = 10
  )

  expect_equal(out$N_pairs[out$Element == "Pb_mgkg"], 20)
  expect_equal(out$Rho[out$Element == "Pb_mgkg"], 1)
  expect_true(is.na(out$Rho[out$Element == "Zn_mgkg"]))
})

test_that("robust_corr_vs_Ag delegates to robust_cor_table", {
  set.seed(1)
  df <- tibble::tibble(
    Ag_mgkg = 1:12,
    Pb_mgkg = 1:12
  )

  out <- robust_corr_vs_Ag(df, metals = "Pb_mgkg", n_boot = 10, min_pairs = 5)

  expect_equal(out$Element, "Pb_mgkg")
  expect_equal(out$Rho, 1)
})

test_that("run_pca_robust returns PCA scores and used columns", {
  df <- tibble::tibble(
    Ag_mgkg = seq(1, 10),
    Pb_mgkg = seq(2, 20, by = 2),
    Zn_mgkg = rep(5, 10)
  )

  out <- run_pca_robust(
    df,
    metal_cols = c("Ag_mgkg", "Pb_mgkg", "Zn_mgkg"),
    min_non_na = 5,
    min_complete = 5
  )

  expect_s3_class(out$pca, "prcomp")
  expect_equal(out$used_cols, c("Ag_mgkg", "Pb_mgkg"))
  expect_equal(nrow(out$scores), 10)
})
