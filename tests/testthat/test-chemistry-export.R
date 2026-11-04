export_fixture <- function() {
  root <- tempfile("geocensor_export_");dir.create(root)
  path <- file.path(root,"chemistry.csv")
  x <- tibble::tibble(lab_id=sprintf("%06d",1:2051),species=rep(c("Ag","Pb"),length.out=2051),
    unit="ppm",qvalue="2.0000",qualifier="",source_data_record="original",
    source_rows="source column",rec_id=as.character(1:2051),note="plain")
  x$qvalue[1001] <- "4.11110";x$rec_id[1001] <- "9007199254740993"
  x$qualifier[c(1501,2001,2031)] <- c("G","L","B")
  x$note[1005] <- "quoted comma,\nsecond line"
  x$species[c(2,4,6)] <- c(" ag ","ag",NA_character_)
  x$qvalue[30] <- NA_character_;x$note[31] <- "NA"
  readr::write_csv(x,path)
  list(root=root,path=path,full=read_usgs_chemistry(path))
}

test_that("disk exports preserve records, positions and verifiable provenance", {
  f <- export_fixture()
  x <- export_usgs_chemistry(f$path,file.path(f$root,"all with spaces"),chunk_size=1000)
  expect_identical(read_usgs_chemistry(x$data_file),f$full)
  expect_identical(readr::read_csv(x$source_rows_file,show_col_types=FALSE)$source_data_record,as.numeric(1:2051))
  expect_equal(x$n_input,2051);expect_equal(x$n_selected,2051)
  expect_equal(x$n_chunks,3);expect_equal(x$max_chunk_rows,1000);expect_null(x$elements)
  info <- readRDS(x$info_file)
  expect_identical(info$columns,names(f$full))
  expect_identical(info$source$md5,unname(tools::md5sum(f$path)))
  expect_equal(info$source$bytes,file.info(f$path)$size)
  expect_identical(info$package_version,as.character(getNamespaceVersion("geoCensor")))
  expect_identical(info$missing_markers,c("","NA"))
  files <- readr::read_csv(x$files_file,show_col_types=FALSE)
  expect_identical(files$file,c("chemistry.csv","source_rows.csv","export_info.rds"))
  expect_identical(files$md5,unname(tools::md5sum(file.path(x$output_dir,files$file))))
  expect_equal(files$bytes,file.info(file.path(x$output_dir,files$file))$size)
  expect_false(dir.exists(file.path(f$root,".all with spaces.geoCensor-lock")))
})

test_that("disk selection is exact and empty selections retain their schema", {
  f <- export_fixture();elements <- c("Pb"," ag ")
  x <- export_usgs_chemistry(f$path,file.path(f$root,"selected"),elements,333)
  rows <- which(f$full$species %in% elements)
  expect_identical(read_usgs_chemistry(x$data_file),f$full[rows,])
  expect_identical(readr::read_csv(x$source_rows_file,show_col_types=FALSE)$source_data_record,as.numeric(rows))
  expect_identical(x$elements,elements)
  absent <- export_usgs_chemistry(f$path,file.path(f$root,"absent"),"not_reported")
  expect_identical(read_usgs_chemistry(absent$data_file),f$full[0,])
  expect_equal(absent$n_selected,0);expect_equal(absent$n_input,2051)
  expect_equal(nrow(readr::read_csv(absent$source_rows_file,show_col_types=FALSE)),0)
  readr::write_csv(f$full[0,],f$path)
  empty <- export_usgs_chemistry(f$path,file.path(f$root,"empty"))
  expect_identical(read_usgs_chemistry(empty$data_file),f$full[0,])
  expect_equal(empty$n_input,0);expect_equal(empty$n_selected,0)
  expect_equal(empty$n_chunks,0);expect_equal(empty$max_chunk_rows,0)
})

test_that("late unselected malformed input leaves no output or reservation", {
  f <- export_fixture();dest <- file.path(f$root,"bad")
  prefix <- c("lab_id,species,qvalue,unit","001,Ag,2,ppm",rep("002,Pb,3,ppm",1000))
  for(last in c("003,Pb,4,ppm,extra","003,Pb,4",'003,Pb,"unterminated,ppm')) {
    writeLines(c(prefix,last),f$path)
    expect_error(export_usgs_chemistry(f$path,dest,"Ag",1000),"Chemistry CSV")
    expect_false(file.exists(dest));expect_false(dir.exists(file.path(f$root,".bad.geoCensor-lock")))
  }
  for(header in c("lab_id,species,qvalue,qvalue,unit","lab_id,species,qvalue,,unit","lab_id,species,qvalue")) {
    writeLines(header,f$path)
    expect_error(export_usgs_chemistry(f$path,dest))
    expect_false(file.exists(dest));expect_false(dir.exists(file.path(f$root,".bad.geoCensor-lock")))
  }
})

test_that("existing outputs and reservations belong to their owners", {
  f <- export_fixture();dest <- file.path(f$root,"existing");dir.create(dest)
  expect_error(export_usgs_chemistry(f$path,dest),"already exists")
  writeLines("keep",file.path(dest,"original.txt"));before <- tools::md5sum(file.path(dest,"original.txt"))
  expect_error(export_usgs_chemistry(f$path,dest),"already exists")
  expect_identical(tools::md5sum(file.path(dest,"original.txt")),before)
  lock <- file.path(f$root,".reserved.geoCensor-lock");dir.create(lock)
  writeLines("another owner",file.path(lock,"owner"));before <- tools::md5sum(file.path(lock,"owner"))
  expect_error(export_usgs_chemistry(f$path,file.path(f$root,"reserved")),"reserved")
  expect_identical(tools::md5sum(file.path(lock,"owner")),before)
  expect_false(file.exists(file.path(f$root,"reserved")))
})

test_that("real source changes before final fingerprint prevent publication", {
  f <- export_fixture();dest <- file.path(f$root,"changed")
  fingerprint <- geoCensor:::.export_source_fingerprint;calls <- 0L
  testthat::local_mocked_bindings(.export_source_fingerprint=function(path) {
    calls <<- calls+1L
    if(calls==2L)cat("extra,Pb,2,ppm,,,,,\n",file=path,append=TRUE)
    fingerprint(path)
  },.package="geoCensor")
  expect_error(export_usgs_chemistry(f$path,dest),"Source file changed")
  expect_false(file.exists(dest));expect_false(dir.exists(file.path(f$root,".changed.geoCensor-lock")))
})

test_that("a destination created by another writer is left intact", {
  f <- export_fixture();dest <- file.path(f$root,"appeared")
  fingerprint <- geoCensor:::.export_source_fingerprint;calls <- 0L
  testthat::local_mocked_bindings(.export_source_fingerprint=function(path) {
    calls <<- calls+1L
    if(calls==2L) {dir.create(dest);writeLines("another writer",file.path(dest,"keep.txt"))}
    fingerprint(path)
  },.package="geoCensor")
  expect_error(export_usgs_chemistry(f$path,dest),"appeared during export")
  expect_identical(readLines(file.path(dest,"keep.txt")),"another writer")
  expect_false(dir.exists(file.path(f$root,".appeared.geoCensor-lock")))
})

test_that("disk export validates explicit arguments before creating output", {
  f <- export_fixture();dest <- file.path(f$root,"invalid")
  for(elements in list(character(),NA_character_,c("Ag","Ag")," ",1)) {
    expect_error(export_usgs_chemistry(f$path,dest,elements),"elements must")
  }
  for(size in list(0,-1,1.5,NA_real_,Inf,numeric(),c(1,2),"2",TRUE,1+1i,2147483648)) {
    expect_error(export_usgs_chemistry(f$path,dest,chunk_size=size),"chunk_size must")
  }
  expect_error(export_usgs_chemistry(dirname(f$path),dest),"local CSV")
  expect_error(export_usgs_chemistry(f$path,file.path(f$root,"missing","child")),"existing parent")
  expect_error(export_usgs_chemistry(f$path,NA_character_),"new directory path")
  expect_false(file.exists(dest))
})
