batch_overlap_input <- function() list(
 chemistry=tibble::tibble(lab_id=rep(sprintf('%03d',1:16),2),
  species=rep(c('Ag','Pb'),each=16),qvalue=rep(as.character(1:16),2),
  unit=rep(c('ppb','unsupported'),each=16),amethod=rep(rep(c('A','B'),each=4),4)),
 metadata=tibble::tibble(lab_id=sprintf('%03d',1:16),state=rep(c('WI',NA_character_),each=8),
  x_m=c(NA,0,10,20,20,30,40,50,100,110,120,130,130,140,150,160),y_m=0))

run_batch_overlap <- function(input=batch_overlap_input(),...) analyze_elements(
 input$chemistry,main_df=input$metadata,group_cols='state',cell_sizes=c(10,100),
 x_col='x_m',y_col='y_m',n_offsets=2,...)

test_that('batch overlap is optional and preserves existing pooled and stratified components', {
 input <- batch_overlap_input();saved <- input
 base <- run_batch_overlap(input)
 audited <- run_batch_overlap(input,audit_method_overlap=TRUE)
 new_names <- c('method_overlap_support','method_coverage','method_overlap')
 expect_named(audited,c(names(base),new_names))
 expect_equal(audited[names(base)],base)
 stratified <- run_batch_overlap(input,stratify_methods=TRUE)
 both <- run_batch_overlap(input,stratify_methods=TRUE,audit_method_overlap=TRUE)
 expect_named(both,c(names(stratified),new_names))
 expect_equal(both[names(stratified)],stratified)
 expect_equal(both[new_names],audited[new_names])
 expect_equal(input,saved)
 expect_true(all(is.finite(base$estimates$geomean[base$estimates$element=='Ag'])))
})

test_that('every batch element/domain audit matches a standalone common-grid reconstruction', {
 input <- batch_overlap_input();result <- run_batch_overlap(input,audit_method_overlap=TRUE)
 expect_equal(nrow(result$method_overlap_support),4)
 expect_equal(nrow(result$method_coverage),16)
 expect_equal(nrow(result$method_overlap),8)
 for (element in c('Ag','Pb')) for (domain in c('WI',NA_character_)) {
  observations <- join_soil_main(build_element(input$chemistry,element),input$metadata,
   extra_cols=c('state','x_m','y_m'))
  select <- if(is.na(domain)) is.na(observations$state) else !is.na(observations$state)&observations$state==domain
  reference <- method_spatial_overlap(observations[select,],c(10,100),'x_m','y_m',n_offsets=2)
  map <- c(method_overlap_support='support',method_coverage='methods',method_overlap='pairs')
  for(component in names(map)) {
   table <- result[[component]]
   match_domain <- if(is.na(domain)) is.na(table$state) else !is.na(table$state)&table$state==domain
   rows <- table[table$element==element&match_domain,setdiff(names(table),c('element','state'))]
   expect_equal(rows,reference[[map[[component]]]])
  }
 }
 expect_true(all(result$method_overlap$n_origins==4))
 expect_true(anyNA(result$method_overlap_support$state))
})

test_that('batch audits retain unavailable chemistry/coordinates as missing rather than observed zero overlap', {
 result <- run_batch_overlap(audit_method_overlap=TRUE)
 pb <- result$method_overlap[result$method_overlap$element=='Pb',]
 expect_true(all(is.na(pb$jaccard_mean)))
 expect_true(all(pb$n_origins_evaluated==0))
 expect_true(all(pb$overlap_status=='unavailable_spatial_data'))
 wi <- result$method_overlap_support[result$method_overlap_support$element=='Ag'&
  !is.na(result$method_overlap_support$state),]
 expect_equal(wi$n_usable_chemistry,8)
 expect_equal(wi$n_spatial,7)
 expect_equal(wi$n_unlocated_usable,1)
 input <- batch_overlap_input();input$metadata$x_m <- NA_real_
 unlocated <- run_batch_overlap(input,audit_method_overlap=TRUE)
 expect_true(all(is.na(unlocated$method_overlap$jaccard_mean)))
 expect_true(all(is.na(unlocated$method_overlap_support$grid_anchor_x)))
 expect_equal(unlocated$method_overlap_support$n_unlocated_usable[
  unlocated$method_overlap_support$element=='Ag'],c(8,8))
})

test_that('unknown method buckets and absent species produce no phantom pairs', {
 input <- batch_overlap_input();input$chemistry$amethod <- NULL
 unknown <- run_batch_overlap(input,audit_method_overlap=TRUE)
 expect_equal(nrow(unknown$method_overlap),0)
 expect_true('method_a'%in%names(unknown$method_overlap))
 expect_true(all(is.na(unknown$method_coverage$method_code)))
 expect_true(all(unknown$method_coverage$method_information=='method_not_supplied'))
 support <- unknown$method_overlap_support
 expect_equal(support$n_unresolved_spatial[support$element=='Ag'&!is.na(support$state)],7)
 expect_equal(support$n_unresolved_spatial[support$element=='Ag'&is.na(support$state)],8)
 absent <- analyze_elements(input$chemistry,elements='Au',main_df=input$metadata,
  cell_sizes=10,x_col='x_m',y_col='y_m',n_offsets=1,audit_method_overlap=TRUE)
 expect_equal(absent$method_overlap_support$n_input,0)
 expect_equal(absent$method_overlap_support$element,'Au')
 expect_true(is.na(absent$method_overlap_support$grid_anchor_x))
 expect_equal(nrow(absent$method_coverage),0)
 expect_equal(nrow(absent$method_overlap),0)
})

test_that('batch overlap validates opt-in, geometry and grouping collisions before returning results', {
 input <- batch_overlap_input()
 for(bad in list(NA,1,'TRUE',logical(),c(TRUE,FALSE)))
  expect_error(run_batch_overlap(input,audit_method_overlap=bad),'audit_method_overlap must be TRUE or FALSE')
 expect_error(analyze_elements(input$chemistry,audit_method_overlap=TRUE),'requires cell_sizes')
 expect_error(analyze_elements(input$chemistry,cell_sizes=10,audit_method_overlap=TRUE),'projected coordinates')
 expect_error(analyze_elements(input$chemistry,main_df=input$metadata,cell_sizes=10,
  x_col='x_m',y_col='x_m',audit_method_overlap=TRUE),'distinct numeric projected')
 input$metadata$grid_anchor_x <- 'domain'
 expect_error(analyze_elements(input$chemistry,main_df=input$metadata,group_cols='grid_anchor_x',
  cell_sizes=10,x_col='x_m',y_col='y_m',audit_method_overlap=TRUE),'Grouping columns conflict with result fields')
})
