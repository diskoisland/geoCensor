test_that("chemistry imports preserve rare qualifiers beyond type-guessing rows", {
  x <- tibble::tibble(lab_id=sprintf("%06d",1:3001),species="Ag",unit="ppm",
    qvalue=rep("2.0000",3001),qualifier=rep("",3001),amethod="AES_HF")
  x$qualifier[c(2001,2501,2801,2901)] <- c("G","L","B","T")
  x$qvalue[2001] <- "10.0000"
  path <- tempfile(fileext=".csv");readr::write_csv(x,path)
  data <- read_usgs_chemistry(path)
  expect_true(all(vapply(data,is.character,logical(1))))
  expect_identical(data$lab_id,x$lab_id)
  expect_identical(data$qvalue,x$qvalue)
  expect_identical(data$qualifier[c(2001,2501,2801,2901)],c("G","L","B","T"))
  expect_true(all(is.na(data$qualifier[-c(2001,2501,2801,2901)])))
  parsed <- build_ag_from_csv(path)
  expect_true(parsed$right_cen[2001]);expect_equal(parsed$right_limit[2001],10000)
  expect_true(is.na(parsed$detected_val[2001]))
  expect_true(parsed$left_cen[2501]);expect_equal(parsed$left_limit[2501],2000)
  expect_identical(parsed$validation_status[c(2801,2901)],rep("unknown_qualifier",2))
})

test_that("numeric text cannot lose precision or acquire an upper suffix", {
  x <- tibble::tibble(lab_id=c("001","002","003","004"),species="Ag",unit="ppm",
    qvalue=c("4.11110","4.1111","10.11111","2.0000000000000001"),
    rec_id=c("9007199254740993","9007199254740994","000003","000004"),
    summary=c("4.11110, AES_HF","4.1111, AES_HF","10.11111, AES_HF","2.0000000000000001, AES_HF"),amethod="AES_HF")
  path <- tempfile(fileext=".csv");readr::write_csv(x,path)
  data <- read_usgs_chemistry(path)
  expect_identical(data,x)
  parsed <- build_ag_from_csv(path)
  expect_identical(parsed$qvalue_raw,x$qvalue)
  expect_identical(parsed$right_cen,c(FALSE,TRUE,TRUE,FALSE))
  expect_equal(parsed$detected_val[1],4111.1)
  expect_equal(parsed$right_limit[2:3],c(4000,10000))
  expect_identical(parsed$rec_id,x$rec_id)
})

test_that("text import declares missing markers and defers value validation", {
  path <- tempfile(fileext=".csv")
  writeLines(c('lab_id,species,qvalue,unit,qualifier,note',
    '001,Ag,,ppm,," keep spaces "','002,Ag,NA,ppm,NA,NA',
    '003,Ag,not_numeric,ppm,,','004,Ag, 2.0000 , ppm ,,0'),path)
  data <- read_usgs_chemistry(path)
  expect_identical(data$qvalue[1:2],rep(NA_character_,2))
  expect_identical(data$note[1]," keep spaces ")
  expect_identical(data$qvalue[4]," 2.0000 ")
  expect_identical(data$unit[4]," ppm ")
  parsed <- build_ag_from_csv(path)
  expect_identical(parsed$validation_status,c("missing_value","missing_value","invalid_value","usable"))
  expect_equal(parsed$detected_val[4],2000)
  expect_equal(nrow(parsed),4)
})

test_that("malformed chemistry CSVs and ambiguous headers stop explicitly", {
  path <- tempfile(fileext=".csv")
  writeLines(c("lab_id,species,qvalue,unit","001,Ag,2,ppm,extra"),path)
  expect_error(suppressWarnings(read_usgs_chemistry(path)),"Chemistry CSV parsing")
  expect_error(suppressWarnings(build_ag_from_csv(path)),"Chemistry CSV parsing")
  writeLines(c("lab_id,species,qvalue,unit","001,Ag,2"),path)
  expect_error(suppressWarnings(read_usgs_chemistry(path)),"Chemistry CSV parsing")
  writeLines(c("lab_id,species,qvalue,qvalue,unit","001,Ag,2,3,ppm"),path)
  expect_error(read_usgs_chemistry(path),"unique")
  writeLines(c("lab_id,species,qvalue,,unit","001,Ag,2,x,ppm"),path)
  expect_error(read_usgs_chemistry(path))
  writeLines(c("lab_id,species,qvalue","001,Ag,2"),path)
  expect_error(read_usgs_chemistry(path),"Missing chemistry columns: unit")
  expect_error(read_usgs_chemistry(dirname(path)),"local CSV")
  expect_error(read_usgs_chemistry(NA_character_),"local CSV")
  expect_error(read_usgs_chemistry(c(path,path)),"local CSV")
})
