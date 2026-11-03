sampling_file_fixture <- function() {
  root <- tempfile("geocensor_sampling_file_"); dir.create(root)
  data <- tibble::tibble(lab_id=c("001"," 002 ","003","004","005","006"),
    state=c("WI","IL","WI","NA","NA","WI"),depth=c(" A ","B","","NA","A","B"),
    x=c("001.0000","1","2.5e0","","NA","Inf"),y=c("2","2","3","4","NaN","5"),
    note=c(" literal NA "," keep spaces ","0","","note","other"))
  input <- file.path(root,"metadata.csv");readr::write_csv(data,input,na="")
  list(root=root,data=data,input=input,output=file.path(root,"audit"))
}

test_that("file workflow preserves raw metadata and explicit coordinate provenance", {
  f <- sampling_file_fixture();before <- tools::md5sum(f$input)
  result <- audit_sampling_file(f$input,f$output,group_cols="state",support_cols="depth",
    x_col="x",y_col="y",include_memberships=TRUE)
  expect_identical(result$input_metadata,f$data)
  expect_identical(result$audit$memberships[names(f$data)],f$data)
  expect_identical(readRDS(file.path(f$output,"audit.rds")),result)
  expect_identical(tools::md5sum(f$input),before)
  expect_equal(nrow(result$audit$coordinate_groups),2)
  expect_equal(result$audit$coordinate_groups$n_analysis_domains[1],2)
  expect_equal(result$audit$memberships$n_records_at_position,c(2L,2L,1L,NA,NA,NA))
  expect_true(all(result$audit$memberships$position_membership_status[4:6]=="unavailable_coordinates"))
  expect_equal(sum(result$audit$fields$n_missing),2)
  conversion <- result$coordinate_conversion
  expect_equal(nrow(conversion),12)
  expect_equal(conversion$raw_value[1],"001.0000")
  expect_equal(conversion$numeric_value[1],1)
  expect_equal(conversion$conversion_status[4:6],c("missing_blank","missing_NA","nonfinite"))
  expect_true(is.nan(conversion$numeric_value[11]))
  expect_setequal(list.files(f$output),c("input_metadata.csv","fields.csv","combinations.csv","locations.csv",
    "memberships.csv","coordinate_groups.csv","coordinate_conversion.csv","audit.rds","run_summary.txt","session_info.txt"))
  expect_match(paste(readLines(file.path(f$output,"run_summary.txt")),collapse="\n"),unname(before),fixed=TRUE)
})

test_that("default and no-coordinate modes export only requested tables", {
  f <- sampling_file_fixture();f$data$x[1] <- "not a coordinate";readr::write_csv(f$data,f$input)
  result <- audit_sampling_file(f$input,f$output,support_cols="depth")
  expect_named(result$audit,c("fields","combinations","locations","available_fields"))
  expect_equal(nrow(result$coordinate_conversion),0)
  expect_false(file.exists(file.path(f$output,"memberships.csv")))
  expect_false(file.exists(file.path(f$output,"coordinate_groups.csv")))
  expect_false(file.exists(file.path(f$output,"locations.csv")))
  expect_false(file.exists(file.path(f$output,"coordinate_conversion.csv")))
  f <- sampling_file_fixture()
  result <- audit_sampling_file(f$input,f$output,support_cols="depth",x_col="x",y_col="y")
  expect_true(file.exists(file.path(f$output,"locations.csv")))
  expect_false(file.exists(file.path(f$output,"memberships.csv")))
  expect_equal(nrow(result$coordinate_conversion),12)
})

test_that("invalid coordinate text fails before creating output", {
  for (value in c("41 N","1,2","0x10","NA?")) {
    f <- sampling_file_fixture();f$data$x[2] <- value;readr::write_csv(f$data,f$input)
    expect_error(audit_sampling_file(f$input,f$output,support_cols="depth",x_col="x",y_col="y"),"nonnumeric.*2")
    expect_false(dir.exists(f$output))
  }
  f <- sampling_file_fixture();f$data$x[1] <- "1e999";readr::write_csv(f$data,f$input)
  result <- audit_sampling_file(f$input,f$output,support_cols="depth",x_col="x",y_col="y",include_memberships=TRUE)
  expect_equal(result$coordinate_conversion$conversion_status[1],"nonfinite")
  expect_true(is.na(result$audit$memberships$recorded_position_group[1]))
})

test_that("nonempty destinations and input files are preserved", {
  f <- sampling_file_fixture()
  audit_sampling_file(f$input,f$output,support_cols="depth")
  files <- list.files(f$output,full.names=TRUE);before <- tools::md5sum(files)
  expect_error(audit_sampling_file(f$input,f$output,support_cols="depth"),"not empty")
  expect_identical(tools::md5sum(files),before)
  expect_error(audit_sampling_file(f$input,f$input,support_cols="depth"),"not a directory")
  empty <- file.path(f$root,"empty");dir.create(empty)
  expect_silent(audit_sampling_file(f$input,empty,support_cols="depth"))
})

test_that("invalid schemas and malformed CSVs never produce partial audits", {
  f <- sampling_file_fixture();f$data$lab_id[2] <- " 001 ";readr::write_csv(f$data,f$input)
  expect_error(audit_sampling_file(f$input,f$output,support_cols="depth"),"nonblank and unique")
  expect_false(dir.exists(f$output))
  f <- sampling_file_fixture();f$data$lab_id[1] <- "NA";readr::write_csv(f$data,f$input)
  expect_error(audit_sampling_file(f$input,f$output,support_cols="depth"),"nonblank and unique")
  expect_false(dir.exists(f$output))
  f <- sampling_file_fixture()
  expect_error(audit_sampling_file(f$input,f$output,support_cols="missing"),"audit columns")
  expect_false(dir.exists(f$output))
  writeLines(c("lab_id,depth","001,A,extra"),f$input)
  expect_error(suppressWarnings(audit_sampling_file(f$input,f$output)),"CSV parsing")
  expect_false(dir.exists(f$output))
})

test_that("custom identifiers and file options are validated without guessing", {
  f <- sampling_file_fixture();names(f$data)[1] <- "sample_key";readr::write_csv(f$data,f$input)
  result <- audit_sampling_file(f$input,f$output,id_col="sample_key",support_cols="depth")
  expect_identical(result$input_metadata$sample_key,f$data$sample_key)
  f <- sampling_file_fixture()
  expect_error(audit_sampling_file(f$input,f$output,include_memberships=TRUE),"requires both coordinate")
  expect_error(audit_sampling_file(f$input,f$output,x_col="x"),"both coordinate")
  expect_error(audit_sampling_file(f$input,f$output,x_col="x",y_col="x"),"distinct columns")
  expect_error(audit_sampling_file(f$input,f$output,include_memberships=NA),"TRUE or FALSE")
  expect_error(audit_sampling_file(f$root,f$output),"local CSV")
  expect_error(audit_sampling_file(f$input,NA_character_),"one directory")
  expect_false(dir.exists(f$output))
})
