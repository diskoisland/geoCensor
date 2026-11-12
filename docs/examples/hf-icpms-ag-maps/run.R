#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(identical(args,'--help') || !length(args)) {
  cat('Usage: Rscript run.R NEW_OUTPUT_DIRECTORY [--refit]\nDefault: redraw verified cells. --refit: recompute every cell from the bundled cohort.\n')
  if(!length(args)) quit(status=1) else quit(status=0)
}
if(length(args)>2L || (length(args)==2L && args[2]!='--refit')) stop('Unexpected arguments; see --help')
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
if(!file.exists(script)) script<-gsub('~+~',' ',script,fixed=TRUE)
here <- dirname(normalizePath(script,mustWork=TRUE));source(file.path(here,'functions.R'))
out <- args[1];refit <- length(args)==2L
if(file.exists(out)) stop('Output already exists; choose a new directory')
if(!dir.exists(dirname(out))) stop('Output parent directory must exist')
input <- hf_read(here);observations <- hf_observations(input$data)
d <- hf_project(observations,input$manifest)
if(nrow(d)!=input$manifest$expected$CONUS_reports) stop('Unexpected mapped report count')
hf_validate_membership(d,input$reference,input$manifest)
if(!dir.create(out)) stop('Cannot create output directory')
# Partial outputs are retained if a later check fails; no inputs are overwritten.
cells <- if(refit) hf_refit(d,input$manifest) else input$reference
checks <- hf_compare(cells,input$reference,input$manifest)
saveRDS(cells,file.path(out,'cell-results.RDS'),compress='xz')
utils::write.csv(d[,c('source_chem_row','medium','ix','iy')],file.path(out,'cell-membership.csv'),row.names=FALSE)
boundaries <- hf_boundaries(input$manifest)
for(medium in names(cells)) {
  utils::write.csv(cells[[medium]],file.path(out,paste0(medium,'-cells.csv')),row.names=FALSE,na='NA')
  for(western in c(FALSE,TRUE)) hf_plot(cells[[medium]],input$manifest,boundaries,medium,western,
    file.path(out,paste0(medium,if(western)'-western.png' else '-conus.png')))
}
jsonlite::write_json(list(status='PASS',mode=if(refit)'fresh_R_refit' else 'redraw_archived_cell_results',
  checks=checks,source_snapshot=input$data$snapshot$source_id,inputs=input$manifest$files,
  geometry=input$manifest[c('crs','cell_width_m','origin_m','circle_radius_m')],
  versions=lapply(c('geoCensor','survival','sf','maps'),function(p)if(requireNamespace(p,quietly=TRUE))
    list(package=p,version=as.character(utils::packageVersion(p))) else list(package=p,version=NA_character_)),
  confidence_intervals=FALSE,interpolation=FALSE),file.path(out,'verification.json'),
  auto_unbox=TRUE,pretty=TRUE,digits=17,na='null')
writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
writeLines('Completed all six maps and source/cell checks.',file.path(out,'COMPLETED.txt'))
message('Completed: ',normalizePath(out))
