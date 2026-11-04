run_analysis_example <- function(script,args) {
  file <- system.file("examples",script,package="geoCensor")
  if(!nzchar(file))file <- testthat::test_path("..","..","inst","examples",script)
  log <- tempfile(fileext=".log")
  previous <- Sys.getenv("R_LIBS",unset=NA_character_)
  on.exit(if(is.na(previous))Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS=previous),add=TRUE)
  Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
  executable <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
  status <- suppressWarnings(system2(executable,
    args=shQuote(c(file,args)),stdout=log,stderr=log))
  list(status=status,log=readLines(log,warn=FALSE))
}

test_that("installed all-element command preserves text groups and unlocated rows", {
  root <- tempfile("geocensor_command_");dir.create(root)
  chem <- file.path(root,"chemistry.csv");metadata <- file.path(root,"metadata.csv")
  readr::write_csv(tibble::tibble(lab_id=sprintf("%03d",1:6),species="Ag",unit="ppm",
    qvalue=c("1","2","3","4.11110","5","6"),qualifier=c("","L","G","","",""),
    rec_id=c("9007199254740993",as.character(2:6)),amethod="MS_HF"),chem)
  readr::write_csv(tibble::tibble(lab_id=sprintf("%03d",1:6),state="001",
    x=c("0","Inf","NA","20","30","40"),y=as.character(1:6)),metadata)
  dest <- file.path(root,"result")
  run <- run_analysis_example("analyze_all_elements.R",c(chem,dest,metadata,"x","y","100","state","--audit-source-precision"))
  expect_equal(run$status,0L,info=paste(run$log,collapse="\n"))
  if(run$status!=0L)return(invisible(NULL))
  x <- readRDS(file.path(dest,"analysis.rds"))
  expect_identical(x$estimates$state,"001")
  expect_equal(x$estimates$n_total,6)
  expect_equal(x$estimates$n_left,1)
  expect_equal(x$estimates$n_right,1)
  expect_equal(x$declustering$n_excluded,2)
  expect_identical(x$source_precision$rec_id[1],"9007199254740993")
  expect_identical(x$source_precision$qvalue[4],"4.11110")
  expect_true(all(is.na(x$declustering$ci_low)))
  # Invalid numeric text must fail before creating a new result directory.
  d <- readr::read_csv(metadata,col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
  d$x[5] <- "40 N";readr::write_csv(d,metadata)
  invalid <- file.path(root,"invalid")
  run <- run_analysis_example("analyze_all_elements.R",c(chem,invalid,metadata,"x","y","100","state"))
  expect_true(run$status!=0L)
  expect_match(paste(run$log,collapse="\n"),"coordinate column.*input rows: 5")
  expect_false(dir.exists(invalid))
})
