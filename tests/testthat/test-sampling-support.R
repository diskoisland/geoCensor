metadata_support <- function() tibble::tibble(
 lab_id=c('a','b','c','d'),state=c('WI','WI','WI','ID'),
 depth=c('0-10 cm','10-20 cm',' ',NA),horizon=c('A','B','A',NA),
 mesh_size=c(' <80 mesh ','<80 mesh','<80 mesh',NA),x=c(1,1,2,NA_real_),y=c(1,1,2,3))

test_that('support audit distinguishes recorded differences from missing support', {
 x<-metadata_support()
 out<-sampling_support_audit(x,group_cols='state',x_col='x',y_col='y')
 depth<-out$fields[out$fields$state=='WI' & out$fields$field=='depth',]
 expect_equal(depth$n_samples,3)
 expect_equal(depth$n_recorded,2)
 expect_equal(depth$n_missing,1)
 expect_equal(depth$n_values,2)
 expect_equal(depth$fraction_recorded,2/3)
 expect_equal(sum(out$combinations$n_samples),4)
 expect_equal(nrow(out$locations),2)
 expect_true(out$locations$mixed_recorded_support[out$locations$x==1])
 expect_false(out$locations$has_support_gaps[out$locations$x==1])
 expect_false(out$locations$mixed_recorded_support[out$locations$x==2])
 expect_true(out$locations$has_support_gaps[out$locations$x==2])
 expect_equal(x,metadata_support())
})

test_that('support audit never calls missing metadata equivalent or merges nearby points', {
 x<-metadata_support()
 x$depth[]<-NA_character_;x$horizon[]<-NA_character_;x$mesh_size[]<-NA_character_
 x$x[2]<-1.00001
 out<-sampling_support_audit(x,x_col='x',y_col='y')
 expect_equal(nrow(out$locations),3)
 expect_true(all(out$locations$has_support_gaps))
 expect_true(all(!out$locations$mixed_recorded_support))
 expect_true(all(out$fields$n_values==0))
})

test_that('support audit validates IDs, fields, and coordinate choices', {
 x<-metadata_support()
 x$lab_id[2]<-' a '
 expect_error(sampling_support_audit(x),'nonblank and unique')
 x<-metadata_support()
 expect_error(sampling_support_audit(x,support_cols='not_present'),'audit columns')
 expect_error(sampling_support_audit(x,x_col='x'),'both coordinate')
 expect_error(sampling_support_audit(x,x_col='x',y_col='x'),'distinct numeric')
 expect_error(sampling_support_audit(x,group_cols='depth',support_cols='depth'),'distinct')
 expect_equal(nrow(sampling_support_audit(x)$locations),0)
 empty<-sampling_support_audit(x[FALSE,],x_col='x',y_col='y')
 expect_equal(nrow(empty$locations),0)
 expect_true(all(empty$fields$n_samples==0))
})
