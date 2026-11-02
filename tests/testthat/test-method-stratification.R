stratified_chemistry <- function() {
  tibble::tibble(lab_id = sprintf('%03d',1:26), species = 'Ag', unit = 'ppb',
    qvalue = c(as.character(1:8),as.character(10*(1:8)),as.character(1:8),'bad','-5'),
    amethod = c(rep(c(' ES_SQ ','es_sq'),4),rep(c('AES_HF','aes_hf '),4),
      rep(c(NA_character_,''),4),'ES_SQ','UNKNOWN'))
}

stratum_row <- function(result, code) {
  x <- result$method_estimates
  x[if (is.na(code)) is.na(x$method_code) else !is.na(x$method_code) & x$method_code == code, ]
}

test_that('method stratification is optional and preserves all pooled results', {
  x <- stratified_chemistry(); saved <- x
  baseline <- analyze_elements(x)
  result <- analyze_elements(x,stratify_methods=TRUE)
  expect_named(baseline,c('estimates','distributions','validation','declustering','joins','methods'))
  expect_named(result,c(names(baseline),'method_estimates'))
  expect_equal(result[names(baseline)],baseline)
  expect_equal(x,saved)
  expect_setequal(result$method_estimates$method_code,c('ES_SQ','AES_HF',NA_character_,'UNKNOWN'))
  expect_equal(sum(result$method_estimates$n_input),26)
  expect_equal(sum(result$method_estimates$n_total),25)
  expect_equal(anyDuplicated(names(result$method_estimates)),0L)
  expect_true(all(result$method_estimates$analysis_scope=='unweighted_method_stratum'))
})

test_that('strata normalize codes and retain independent information-review counts', {
  result <- analyze_elements(stratified_chemistry(),stratify_methods=TRUE)
  semi <- stratum_row(result,'ES_SQ'); other <- stratum_row(result,'AES_HF')
  expect_equal(semi$n_input,9)
  expect_equal(semi$n_excluded,1)
  expect_equal(semi$n_total,8)
  expect_equal(semi$geomean,exp(mean(log(1:8))),tolerance=1e-6)
  expect_equal(other$geomean,10*semi$geomean,tolerance=1e-6)
  expect_true(semi$review_semiquantitative)
  expect_equal(semi$method_information,'review_semiquantitative')
  expect_match(semi$method_review_source,'pubs.usgs.gov',fixed=TRUE)
  expect_equal(other$method_information,'method_not_reviewed')
  expect_true(is.na(other$review_semiquantitative))
  expect_true(is.na(other$method_review_source))
  expect_equal(semi$review_status,'review_information')
  cleared <- analyze_elements(stratified_chemistry(),stratify_methods=TRUE,min_exact_review=8)
  expect_equal(stratum_row(cleared,'ES_SQ')$review_status,'counts_above_review_thresholds')
  expect_true(stratum_row(cleared,'ES_SQ')$review_semiquantitative)
})

test_that('missing and unspecified buckets remain explicitly unresolved', {
  result <- analyze_elements(stratified_chemistry(),stratify_methods=TRUE)
  missing <- stratum_row(result,NA_character_)
  unknown <- stratum_row(result,'UNKNOWN')
  expect_equal(missing$n_input,8)
  expect_equal(missing$method_information,'method_not_supplied')
  expect_equal(missing$n_method_unspecified_used,8)
  expect_equal(missing$method_review_status,'method_information_unresolved')
  expect_true(is.na(missing$review_semiquantitative))
  expect_equal(unknown$n_input,1)
  expect_equal(unknown$method_information,'method_unspecified')
  expect_equal(unknown$method,'insufficient_data')
  expect_true(is.na(unknown$geomean))
  expect_equal(unknown$n_left,1)
  x <- stratified_chemistry(); x$amethod <- NULL
  absent <- analyze_elements(x,stratify_methods=TRUE)$method_estimates
  expect_equal(nrow(absent),1)
  expect_true(is.na(absent$method_code))
  expect_equal(absent$n_method_unspecified_used,25)
})

test_that('stratification stays within element and domain and uses chemistry rather than coordinates', {
  x <- stratified_chemistry(); other <- x
  other$lab_id <- paste0(other$lab_id,'b'); other$qvalue <- rep('10',nrow(other))
  chemistry <- dplyr::bind_rows(x,other)
  metadata <- tibble::tibble(lab_id=chemistry$lab_id,
    domain=rep(c('first','second'),each=nrow(x)),
    x_m=c(NA,rep(0,nrow(chemistry)-1L)),y_m=0)
  more <- chemistry; more$species <- 'Pb'; more$qvalue <- '20'
  chemistry <- dplyr::bind_rows(chemistry,more)
  result <- analyze_elements(chemistry,main_df=metadata,group_cols='domain',
    cell_sizes=20,x_col='x_m',y_col='y_m',n_offsets=2,stratify_methods=TRUE)
  expect_equal(nrow(result$method_estimates),16)
  rows <- result$method_estimates
  first <- rows[rows$element=='Ag' & rows$domain=='first' & !is.na(rows$method_code) & rows$method_code=='ES_SQ',]
  expect_equal(first$n_total,8)
  expect_true(all(is.na(result$declustering$ci_low)))
  for (element in c('Ag','Pb')) for (domain in c('first','second')) {
    selected <- rows$element==element & rows$domain==domain
    ids <- metadata$lab_id[metadata$domain==domain]
    ref <- analyze_elements(chemistry[chemistry$lab_id%in%ids,],elements=element,stratify_methods=TRUE)$method_estimates
    expect_equal(rows[selected,setdiff(names(rows),'domain')],ref)
  }
})

test_that('unusable and absent strata never create fitted values or phantom methods', {
  x <- stratified_chemistry(); x$unit[x$amethod=='ES_SQ' & !is.na(x$amethod)] <- 'unsupported'
  result <- analyze_elements(x,elements=c('Ag','Au'),stratify_methods=TRUE)
  expect_false(any(result$method_estimates$element=='Au'))
  empty <- analyze_elements(x,elements='Au',stratify_methods=TRUE)$method_estimates
  expect_equal(nrow(empty),0)
  expect_true(all(c('method_code','geomean','method_information')%in%names(empty)))
  x <- stratified_chemistry(); x$unit <- 'unsupported'
  unusable <- analyze_elements(x,stratify_methods=TRUE)$method_estimates
  expect_equal(nrow(unusable),4)
  expect_true(all(is.na(unusable$geomean)))
  expect_true(all(unusable$n_total==0))
  expect_true(all(unusable$method_review_status=='no_usable_data'))
})

test_that('failed stratum fits stay missing and method concerns remain visible', {
  original <- survival::survreg
  testthat::local_mocked_bindings(survreg=function(...) {
    warning('did not converge'); original(...)
  },.package='survival')
  result <- analyze_elements(stratified_chemistry(),stratify_methods=TRUE)
  semi <- stratum_row(result,'ES_SQ')
  expect_equal(semi$method,'failed_fit')
  expect_true(is.na(semi$geomean))
  expect_true(semi$review_semiquantitative)
  expect_equal(semi$n_total,8)
})

test_that('method stratification validates its option and prevents output-name collisions', {
  x <- stratified_chemistry()
  for (bad in list(NA,1,'TRUE',logical(),c(TRUE,FALSE)))
    expect_error(analyze_elements(x,stratify_methods=bad),'stratify_methods must be TRUE or FALSE')
  x$method_information <- 'recorded domain'
  expect_error(analyze_elements(x,group_cols='method_information',stratify_methods=TRUE),
    'Grouping columns conflict with result fields')
})
