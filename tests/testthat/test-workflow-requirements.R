test_that("workflow guards accept final versions and reject missing candidate capabilities", {
  helper <- system.file("examples", "workflow_requirements.R", package = "geoCensor")
  if (!nzchar(helper)) helper <- testthat::test_path("..", "..", "inst", "examples", "workflow_requirements.R")
  env <- new.env(parent = baseenv())
  sys.source(helper, envir = env)
  for (workflow in c("all", "single", "audit", "export", "illinois")) {
    expect_true(env$geocensor_require_workflow(workflow, version = "1.2.0"))
    expect_true(env$geocensor_require_workflow(workflow, version = "1.2.0.9008"))
    expect_error(env$geocensor_require_workflow(workflow, version = "1.1.0"), "final geoCensor 1.2.0")
  }
  stale <- new.env(parent = emptyenv())
  stale$export_usgs_chemistry <- function(chem_csv, output_dir) NULL
  expect_error(env$geocensor_require_workflow("export", version = "1.2.0", namespace = stale),
    "Missing workflow capabilities: export_usgs_chemistry")
  stale$export_usgs_chemistry <- function(chem_csv, output_dir, elements, chunk_size) NULL
  expect_true(env$geocensor_require_workflow("export", version = "1.2.0", namespace = stale))
  expect_error(env$geocensor_require_workflow("all", version = "1.2.0", namespace = stale),
    "Missing workflow capabilities")
})
