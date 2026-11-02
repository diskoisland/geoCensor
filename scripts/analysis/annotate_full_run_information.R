#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr);library(dplyr);pkgload::load_all('.',quiet=TRUE)})
source<-'output/usgs_all_elements_2026-10-31'
folder<-'output/usgs_information_review_2026-10-31'
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
summary<-list()
for(medium in c('soil','sed')){
 for(scope in c('pooled','state','spatial')){
  filename<-switch(scope,pooled='pooled_estimates.csv',state='state_estimates.csv',spatial='state_declustering.csv')
  x<-read_csv(file.path(source,medium,filename),show_col_types=FALSE)
  out<-flag_estimate_information(x,min_exact=10,min_exact_fraction=.05)
  stopifnot(identical(out$geomean,x$geomean),nrow(out)==nrow(x))
  write_csv(out,file.path(folder,paste0(medium,'_',scope,'_information.csv')))
  summary[[paste(medium,scope)]]<-out%>%count(.data$review_status,name='n_estimates')%>%
   mutate(medium=medium,scope=scope)
 }
}
counts<-bind_rows(summary)
write_csv(counts,file.path(folder,'review_counts.csv'))
writeLines(c(
 '# Count-based information review of the completed USGS run','',
 'These files annotate existing results; they do not recalculate, filter, or replace estimates.',
 'The default review conventions flag fewer than 10 exact observations or less than 5% exact observations.',
 'These are configurable review thresholds, not validated requirements for reliable inference.',
 'Exact counts exclude both left and right censoring; a right-censored value is not classified as a nondetection.',
 'Missing estimates remain estimate_unavailable. Passing the count flags does not establish model fit, spatial independence, support comparability, or freedom from preferential selection.',
 'Each spatial row is one element/domain/cell-size combination, so counts across spatial rows must not be interpreted as counts of independent analyses.',
 'Original source and full-run result files are preserved. All annotated geometric means were checked to equal the originals.',
 'Reproduce with scripts/analysis/annotate_full_run_information.R.'
),file.path(folder,'README.md'))
print(counts)
