test_that("selection preserves fields, order and parsed record ordinals", {
  x <- tibble::tibble(lab_id=sprintf("%06d",1:2051),
    species=rep(c("Ag","Pb"),length.out=2051),unit="ppm",qvalue="2.0000",
    qualifier="",rec_id=as.character(1:2051),source_rows="original column",
    source_record_number="000007",note="plain")
  x$rec_id[1001] <- "9007199254740993"; x$qvalue[1001] <- "4.11110"
  x$qualifier[c(1501,2001,2031)] <- c("G","L","B")
  x$note[1005] <- "quoted comma,\nsecond line"
  path <- tempfile(fileext=".csv"); readr::write_csv(x,path)
  full <- read_usgs_chemistry(path)
  selected <- read_usgs_chemistry_select(path,"Ag",1000)
  rows <- which(full$species=="Ag")
  expect_identical(selected$data,full[rows,])
  expect_identical(selected$source_rows,as.numeric(rows))
  expect_equal(selected$n_input,2051); expect_equal(selected$n_chunks,3)
  expect_equal(selected$max_chunk_rows,1000)
  expect_identical(selected$data$note[selected$source_rows==1005],"quoted comma,\nsecond line")
  expect_identical(selected$data$rec_id[selected$source_rows==1001],"9007199254740993")
  expect_identical(selected$data$qvalue[selected$source_rows==1001],"4.11110")
  expect_identical(selected$data$source_rows,rep("original column",length(rows)))
  parsed <- build_element(selected$data,"Ag")
  expect_true(parsed$right_cen[selected$source_rows==1501])
  expect_true(parsed$left_cen[selected$source_rows==2001])
  expect_identical(parsed$validation_status[selected$source_rows==2031],"unknown_qualifier")
  multiple <- read_usgs_chemistry_select(path,c("Pb","Ag"),333)
  expect_identical(multiple$data,full)
  expect_identical(multiple$source_rows,as.numeric(seq_len(nrow(full))))
  expect_equal(multiple$max_chunk_rows,333)
  expect_identical(multiple$elements,c("Pb","Ag"))
})

test_that("matching is exact and empty selections preserve the schema", {
  path <- tempfile(fileext=".csv")
  writeLines(c("lab_id,species,qvalue,unit,note","001,Ag,NA,ppm,NA",
    '002, ag , 2.0000 , ppm ," keep "',"",'003,ag,3,ppm,"NA"',
    "004,,4,ppm,"),path)
  full <- read_usgs_chemistry(path)
  selected <- read_usgs_chemistry_select(path,c("ag"," ag "),1)
  expect_identical(selected$data,full[c(2,3),])
  expect_identical(selected$source_rows,c(2,3)); expect_equal(selected$n_input,4)
  expect_identical(selected$data$qvalue[1]," 2.0000 ")
  expect_identical(selected$data$note[2],NA_character_)
  absent <- read_usgs_chemistry_select(path,"Cu",2)
  expect_identical(absent$data,full[0,]); expect_identical(absent$source_rows,numeric())
  expect_equal(absent$n_input,4)
  writeLines("lab_id,species,qvalue,unit,note",path)
  empty <- read_usgs_chemistry_select(path,"Ag")
  expect_identical(empty$data,read_usgs_chemistry(path))
  expect_equal(empty$n_input,0); expect_equal(empty$n_chunks,0)
  expect_equal(empty$max_chunk_rows,0)
})

test_that("late malformed unselected records cannot yield partial output", {
  path <- tempfile(fileext=".csv")
  prefix <- c("lab_id,species,qvalue,unit","001,Ag,2,ppm",rep("002,Pb,3,ppm",1000))
  for(last in c("003,Pb,4,ppm,extra","003,Pb,4",'003,Pb,"unterminated,ppm')) {
    writeLines(c(prefix,last),path)
    expect_error(read_usgs_chemistry_select(path,"Ag",1000),"Chemistry CSV parsing")
  }
  writeLines(c("lab_id,species,qvalue,qvalue,unit","001,Ag,2,3,ppm"),path)
  expect_error(read_usgs_chemistry_select(path,"Ag"),"unique")
  writeLines(c("lab_id,species,qvalue,,unit","001,Ag,2,x,ppm"),path)
  expect_error(read_usgs_chemistry_select(path,"Ag"))
  writeLines(c("lab_id,species,qvalue","001,Ag,2"),path)
  expect_error(read_usgs_chemistry_select(path,"Ag"),"Missing chemistry columns: unit")
})

test_that("selective reader validates path and explicit options", {
  path <- tempfile(fileext=".csv"); writeLines("lab_id,species,qvalue,unit",path)
  for(value in list(character(),NA_character_,c("Ag","Ag"),c("Ag",NA)," ",1)) {
    expect_error(read_usgs_chemistry_select(path,value),"elements must")
  }
  for(value in list(0,-1,1.5,NA_real_,Inf,numeric(),c(1,2),"2",TRUE,2147483648)) {
    expect_error(read_usgs_chemistry_select(path,"Ag",value),"chunk_size must")
  }
  for(value in list(NA_character_,character(),c(path,path),"",dirname(path),paste0(path,"missing"))) {
    expect_error(read_usgs_chemistry_select(value,"Ag"),"local CSV")
  }
})
