method_example <- function() {
 tibble::tibble(lab_id=sprintf('%03d',1:10),species='Ag',unit='ppm',
  qvalue=c('10','20','30','40','50','60','-5','100.1111','bad','0'),
  amethod=c(' ES_SQ ','es_h2o_p','AES_HF',NA,'','UNKNOWN','ES_SQ','ES_SQ','ES_SQ','ES_SQ'))
}

test_that('method review preserves values and unknown methods stay unresolved', {
 x<-tibble::tibble(amethod=factor(c(' ES_SQ ','es_h2o_p','AES_HF','UNKNOWN','',NA,'custom')),
                   qvalue=c('1','2','3','4','-1','5.1111','6'))
 out<-flag_analytical_methods(x)
 expect_equal(out[names(x)],x)
 expect_equal(out$review_semiquantitative,c(TRUE,TRUE,rep(NA,5)))
 expect_equal(out$method_review_status,c('review_semiquantitative','review_semiquantitative',
  'method_not_reviewed','method_unspecified','method_not_supplied','method_not_supplied','method_not_reviewed'))
 expect_true(all(grepl('pubs.usgs.gov',out$method_review_source[1:2],fixed=TRUE)))
 expect_true(all(is.na(out$method_review_source[3:7])))
 expect_equal(flag_analytical_methods(out),out)
 expect_true(all(is.na(flag_analytical_methods(data.frame(x=1:2))$review_semiquantitative)))
 empty<-flag_analytical_methods(data.frame(amethod=character()))
 expect_equal(nrow(empty),0)
 expect_type(empty$review_semiquantitative,'logical')
 expect_type(empty$method_review_status,'character')
})

test_that('custom method columns are supported without guessing invalid inputs', {
 x<-data.frame(assay='ES_SQ',value=10)
 expect_true(flag_analytical_methods(x,'assay')$review_semiquantitative)
 expect_equal(flag_analytical_methods(x,'assay')$assay,x$assay)
 for(bad in list(NA_character_,character(),c('a','b'),'','review_semiquantitative'))
  expect_error(flag_analytical_methods(x,bad),'method_col')
 expect_error(flag_analytical_methods(1:3),'data frame')
 expect_error(flag_analytical_methods(data.frame(amethod=I(list('ES_SQ')))),'atomic vector')
})

test_that('method concerns are separate from chemistry validity and coding', {
 x<-method_example()
 out<-build_element(x,'Ag')
 expect_equal(out$lab_id,x$lab_id)
 expect_equal(out$validation_status,c(rep('usable',8),'invalid_value','nonpositive_detection'))
 expect_equal(out$detected_val[1],10000)
 expect_equal(out$left_limit[7],5000)
 expect_equal(out$right_limit[8],100000)
 expect_true(out$review_semiquantitative[9])
 expect_true(is.na(out$detected_val[9]))
 expect_false(any(c('interval_lower','interval_upper')%in%names(out)))
})

test_that('batch counts distinguish input flags from records actually used', {
 x<-method_example()
 out<-analyze_elements(x)
 e<-out$estimates
 expect_equal(e$n_total,8)
 expect_equal(e$n_exact,6)
 expect_equal(e$n_semiquantitative_input,6)
 expect_equal(e$n_semiquantitative_used,4)
 expect_equal(e$n_method_unspecified_used,3)
 expect_equal(e$n_method_not_reviewed_used,1)
 expect_true(e$review_semiquantitative)
 expect_equal(e$method_review_status,'review_semiquantitative')
 expect_equal(anyDuplicated(names(e)),0L)
 expect_equal(tail(names(e),6),c('n_semiquantitative_input','n_semiquantitative_used',
  'n_method_unspecified_used','n_method_not_reviewed_used','review_semiquantitative','method_review_status'))
 expect_equal(sum(out$methods$n_input),10)
 expect_equal(sum(out$methods$n_usable),8)
 expect_equal(sum(out$methods$n_uncensored),6)
 expect_equal(sum(out$methods$n_left),1)
 expect_equal(sum(out$methods$n_right),1)
 expect_true(all(out$methods$analysis_scope=='unweighted_estimate'))
 reference<-fit_cen_lognorm(build_element(x,'Ag'))
 expect_equal(e[c('geomean','ci_low','ci_high')],reference[c('geomean','ci_low','ci_high')])
})

test_that('unknown, excluded and absent data never create a clear method flag', {
 x<-method_example()
 x$amethod<-'AES_HF'
 e<-analyze_elements(x)$estimates
 expect_true(is.na(e$review_semiquantitative))
 expect_equal(e$method_review_status,'method_information_unresolved')
 absent<-analyze_elements(x,elements='Au')
 expect_equal(absent$estimates$method_review_status,'no_usable_data')
 expect_true(is.na(absent$estimates$review_semiquantitative))
 expect_equal(nrow(absent$methods),0)
 bad<-x[9:10,];bad$amethod<-'ES_SQ'
 e<-analyze_elements(bad)$estimates
 expect_equal(e$n_semiquantitative_input,2)
 expect_equal(e$n_semiquantitative_used,0)
 expect_true(is.na(e$review_semiquantitative))
 expect_equal(e$method_review_status,'no_usable_data')
})

test_that('method review does not disappear when count review thresholds clear', {
 x<-tibble::tibble(lab_id=sprintf('%03d',1:20),species='Ag',qvalue=as.character(1:20),unit='ppm',amethod='ES_SQ')
 e<-analyze_elements(x)$estimates
 expect_equal(e$review_status,'counts_above_review_thresholds')
 expect_true(e$review_semiquantitative)
 expect_equal(e$n_semiquantitative_used,20)
})

test_that('spatial method counts use the coordinate subset instead of all chemistry', {
 x<-method_example()[1:8,]
 x$amethod<-c('ES_SQ',rep('AES_HF',7))
 metadata<-tibble::tibble(lab_id=x$lab_id,state='WI',x_m=c(NA,0,0,0,10,20,30,40),y_m=0)
 result<-analyze_elements(x,main_df=metadata,group_cols='state',cell_sizes=c(5,10),x_col='x_m',y_col='y_m')
 expect_equal(result$estimates$n_semiquantitative_used,1)
 expect_true(result$estimates$review_semiquantitative)
 expect_equal(result$declustering$n_semiquantitative_input,c(1,1))
 expect_equal(result$declustering$n_semiquantitative_used,c(0,0))
 expect_true(all(is.na(result$declustering$review_semiquantitative)))
 expect_equal(result$declustering$n_total,c(7,7))
 expect_true(all(is.na(result$declustering$ci_low)))
 expect_equal(anyDuplicated(names(result$declustering)),0L)
})


test_that('grouping by raw methods works and colliding output names are rejected', {
 x<-method_example()
 out<-analyze_elements(x,group_cols='amethod')
 expect_equal(sum(out$estimates$n_input),nrow(x))
 expect_true('amethod'%in%names(out$methods))
 expect_error(analyze_elements(x,group_cols='method_review_status'),
  'Grouping columns conflict with result fields')
})
