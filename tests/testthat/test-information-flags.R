test_that('information flags preserve estimates and distinguish exact data from censoring', {
 x<-tibble::tibble(n_total=c(100,100,100,0),n_left=c(95,70,30,0),n_right=c(0,20,20,0),geomean=c(2,3,4,NA_real_))
 out<-flag_estimate_information(x)
 expect_equal(out$geomean,x$geomean)
 expect_equal(out$n_exact,c(5,10,50,0))
 expect_equal(out$exact_fraction,c(.05,.1,.5,NA_real_))
 expect_equal(out$review_few_exact,c(TRUE,FALSE,FALSE,TRUE))
 expect_equal(out$review_status,c('review_information','counts_above_review_thresholds','counts_above_review_thresholds','estimate_unavailable'))
 expect_true(is.na(out$review_low_exact_fraction[4]))
 relaxed<-flag_estimate_information(x,min_exact=0,min_exact_fraction=0)
 expect_equal(relaxed$review_status[1],'counts_above_review_thresholds')
 expect_equal(flag_estimate_information(out),out)
})

test_that('information flags reject inconsistent counts and invalid thresholds', {
 x<-tibble::tibble(n_total=10,n_left=5,n_right=6,geomean=1)
 expect_error(flag_estimate_information(x),'exceed')
 x$n_right<-1
 expect_error(flag_estimate_information(x,min_exact=NA_real_),'nonnegative')
 expect_error(flag_estimate_information(x,min_exact_fraction=1.1),'\\[0, 1\\]')
 x$n_left<-NA_real_
 expect_error(flag_estimate_information(x),'finite nonnegative')
})

test_that('batch analyses add review flags without filtering sparse estimates', {
 x<-tibble::tibble(lab_id=as.character(1:8),species='Ag',qvalue=as.character(1:8),unit='ppm')
 out<-analyze_elements(x)$estimates
 expect_equal(out$n_exact,8)
 expect_equal(out$review_status,'review_information')
 expect_true(is.finite(out$geomean))
 out2<-analyze_elements(x,min_exact_review=5)$estimates
 expect_equal(out2$geomean,out$geomean)
 expect_equal(out2$review_status,'counts_above_review_thresholds')
})
