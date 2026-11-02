#!/usr/bin/env Rscript
# Fresh observational validation of the checked method-aware workflow.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop('Usage: Rscript run_usgs_method_validation.R output_dir pilot|expansion soil|sed')
out <- normalizePath(args[1],mustWork=TRUE);stage <- args[2];medium <- args[3]
stopifnot(stage%in%c('pilot','expansion'),medium%in%c('soil','sed'),stage!='pilot'||medium=='soil',
 file.exists(file.path(out,'PREPARED.txt')))
suppressPackageStartupMessages({library(data.table);library(sf)})
setDTthreads(2L);Sys.setenv(TZ='America/Chicago')
config <- read.csv(file.path(out,'config.csv'),colClasses='character',check.names=FALSE)
value <- function(key)config$value[match(key,config$key)]
.libPaths(c(value('library_dir'),.libPaths()))
suppressPackageStartupMessages(library(geoCensor))
stopifnot(as.character(packageVersion('geoCensor'))==value('package_version'))
widths <- as.numeric(strsplit(value('cell_sizes_m'),',',fixed=TRUE)[[1]])
n_offsets <- as.integer(value('n_offsets'))
groups <- strsplit(value('group_cols'),',',fixed=TRUE)[[1]]
source_manifest <- read.csv(file.path(out,'source_manifest.csv'),colClasses='character')
input_manifest <- read.csv(file.path(out,'input_manifest.csv'),colClasses='character')
verify_manifest <- function(m) {
 actual <- unname(tools::md5sum(m$file))
 if(anyNA(actual)||!identical(actual,m$md5))stop('Frozen source/input changed; use a fresh validation directory')
}
verify_manifest(source_manifest);verify_manifest(input_manifest)
prepared <- file.path(out,'inputs',medium)
raw <- as.data.frame(fread(file.path(prepared,'main_selected.csv'),colClasses='character',na.strings=NULL))
chemistry <- as.data.frame(fread(file.path(prepared,'chemistry_selected.csv'),colClasses='character',na.strings=NULL))
stopifnot(all(c('amethod','summary','parameter','mthd_parameter','source_record_number')%in%names(chemistry)))
raw$lab_id <- trimws(raw$lab_id)
stopifnot(!anyNA(raw$lab_id),all(nzchar(raw$lab_id)),!anyDuplicated(raw$lab_id),
 all(toupper(trimws(raw$state))==value('state')),
 all(toupper(trimws(raw$country))%in%c('US','USA','UNITED STATES','U.S.')))
old <- file.path('output','usgs_all_elements_2026-10-31',medium)
projected <- as.data.frame(fread(file.path(old,'metadata_projected.csv'),colClasses=c(lab_id='character')))
projected$lab_id <- trimws(projected$lab_id)
stopifnot(!anyDuplicated(projected$lab_id))
row_index <- match(raw$lab_id,projected$lab_id);stopifnot(!anyNA(row_index))
p <- projected[row_index,,drop=FALSE]
stopifnot(all(toupper(trimws(p$state))==value('state')),
 all(toupper(trimws(p$country))%in%c('US','USA','UNITED STATES','U.S.')))
lat <- suppressWarnings(as.numeric(raw$lat_wgs84));lon <- suppressWarnings(as.numeric(raw$long_wgs84))
valid <- is.finite(lat)&abs(lat)<=90&is.finite(lon)&abs(lon)<=180
stopifnot(identical(valid,is.finite(p$x_m)&is.finite(p$y_m)))
if(any(valid))stopifnot(max(abs(lat[valid]-p$lat_wgs84[valid]))<1e-10,
 max(abs(lon[valid]-p$long_wgs84[valid]))<1e-10)
definitions <- as.data.frame(fread(file.path(old,'projection_definitions.csv')))
definition <- definitions[definitions$state==value('state'),,drop=FALSE]
stopifnot(nrow(definition)==1L,definition$n_coordinates==sum(valid))
max_error <- NA_real_
if(any(valid)) {
 points <- st_as_sf(data.frame(lon=lon[valid],lat=lat[valid]),coords=c('lon','lat'),crs=4326)
 reproduced <- st_coordinates(st_transform(points,definition$crs[1]))
 max_error <- max(abs(reproduced[,1]-p$x_m[valid]),abs(reproduced[,2]-p$y_m[valid]))
 stopifnot(max_error<=1e-5)
}
metadata <- raw;metadata$state <- toupper(trimws(metadata$state))
metadata$method_collected <- trimws(metadata$method_collected)
metadata$method_collected[metadata$method_collected==''] <- NA_character_
metadata$x_m <- p$x_m;metadata$y_m <- p$y_m
context <- file.path(out,'context',medium)
context_signature <- list(source_md5=unname(tools::md5sum(source_manifest$file)),
 input_md5=unname(tools::md5sum(input_manifest$file)),
 protocol=unname(tools::md5sum(file.path(out,'PROTOCOL.md'))),
 R=R.version.string,sf=as.character(packageVersion('sf')),sf_external=sf_extSoftVersion())
if(file.exists(file.path(context,'COMPLETED.txt'))) {
 stopifnot(identical(readRDS(file.path(context,'signature.rds')),context_signature),
  identical(readRDS(file.path(context,'metadata.rds')),metadata))
} else {
 if(dir.exists(context)&&length(list.files(context,all.files=TRUE,no..=TRUE)))stop('Incomplete context exists; preserve and use a fresh run')
 dir.create(context,recursive=TRUE)
 saveRDS(context_signature,file.path(context,'signature.rds'));saveRDS(metadata,file.path(context,'metadata.rds'))
 fwrite(metadata,file.path(context,'metadata.csv'),na='NA')
 fwrite(definition,file.path(context,'projection_definition.csv'))
 fwrite(data.frame(medium=medium,n_metadata=nrow(raw),n_finite_coordinates=sum(valid),
  max_reprojection_error_m=max_error,tolerance_m=1e-5),file.path(context,'projection_reconciliation.csv'))
 support_fields <- intersect(c('project_name','primary_class','depth','mesh_size','prep','collection_protocol',
  'sample_source','geologic_age','stratigraphy','mineralization','date_collect','coord_prec_est',
  'previous_lab_id','sample_zone','horizon'),names(metadata))
 audit <- sampling_support_audit(metadata,group_cols=groups,support_cols=support_fields,x_col='x_m',y_col='y_m')
 for(name in c('fields','combinations','locations'))fwrite(audit[[name]],file.path(context,paste0('support_',name,'.csv')),na='NA')
 coverage <- lapply(c(support_fields,'method_collected'),function(field) {
  values <- trimws(as.character(metadata[[field]]));missing <- is.na(values)|values==''
  data.frame(medium=medium,field=field,n_metadata=nrow(metadata),n_missing=sum(missing),
   n_literal_unknown=sum(!missing&tolower(values)=='unknown'),n_other_recorded=sum(!missing&tolower(values)!='unknown'))
 })
 fwrite(rbindlist(coverage),file.path(context,'field_coverage.csv'))
 fwrite(data.frame(medium=medium,reporting_design_status='unestablished',
  geological_population_status='unestablished',physical_site_status='unestablished',
  spatial_interval_status='experimental',p_value=NA_real_),file.path(context,'inference_status.csv'),na='NA')
 writeLines('Metadata/projection reconciliation complete; labels do not establish support equivalence.',file.path(context,'COMPLETED.txt'))
}
if(stage=='pilot') {
 metadata <- metadata[!is.na(metadata$method_collected)&metadata$method_collected==value('pilot_collection'),,drop=FALSE]
 elements <- strsplit(value('pilot_elements'),',',fixed=TRUE)[[1]]
} else {
 elements <- sort(unique(chemistry$species))
 if(anyNA(elements)||any(!nzchar(trimws(elements))))stop('Missing species labels require review; do not invent an element')
}
chemistry <- chemistry[trimws(chemistry$lab_id)%in%metadata$lab_id,,drop=FALSE]
if(anyDuplicated(chemistry[c('lab_id','species')]))stop('Duplicate sample/species keys; no averaging permitted')
folder <- file.path(out,'results',stage,medium);signature_file <- file.path(folder,'signature.rds')
installed <- system.file(package='geoCensor')
signature <- list(stage=stage,medium=medium,elements=elements,cell_sizes=widths,n_offsets=n_offsets,groups=groups,
 n_metadata=nrow(metadata),selected_inputs=tools::md5sum(file.path(prepared,c('main_selected.csv','chemistry_selected.csv'))),
 context=context_signature,config=unname(tools::md5sum(file.path(out,'config.csv'))),
 installed=tools::md5sum(file.path(installed,c('DESCRIPTION','NAMESPACE','R/geoCensor.rdb','R/geoCensor.rdx'))),
 R=R.version.string,packages=vapply(c('geoCensor','data.table','sf','survival','dplyr','readr'),
  function(package)as.character(packageVersion(package)),character(1)))
if(file.exists(signature_file)) {
 stopifnot(identical(readRDS(signature_file),signature))
 if(file.exists(file.path(folder,'COMPLETED.txt'))) {
  verify_manifest(source_manifest);verify_manifest(input_manifest)
  cat(stage,medium,'completed matching run; unchanged\n');quit(status=0L)
 }
} else {
 if(dir.exists(folder)&&length(list.files(folder,all.files=TRUE,no..=TRUE)))stop('Use a fresh result directory')
 dir.create(file.path(folder,'checkpoints'),recursive=TRUE);saveRDS(signature,signature_file)
 writeLines(capture.output(sessionInfo()),file.path(folder,'sessionInfo.txt'))
 fwrite(data.frame(stage=stage,medium=medium,n_selected_metadata=nrow(metadata),
  n_selected_chemistry=nrow(chemistry),n_requested_species=length(elements)),file.path(folder,'cohort.csv'))
}
for(i in seq_along(elements)) {
 element <- elements[i];checkpoint <- file.path(folder,'checkpoints',sprintf('%03d.rds',i))
 if(file.exists(checkpoint)) {
  cached <- readRDS(checkpoint);stopifnot(identical(cached$element,element),identical(cached$signature,signature))
  cat(stage,medium,i,'/',length(elements),element,'cached\n');flush.console();next
 }
 d <- chemistry[chemistry$species==element,,drop=FALSE]
 d$lab_id_raw <- d$lab_id
 cat(stage,medium,i,'/',length(elements),element,'records',nrow(d),'starting\n');flush.console()
 start <- proc.time()[3];warnings <- character()
 result <- withCallingHandlers(analyze_elements(d,elements=element,main_df=metadata,group_cols=groups,
  cell_sizes=widths,x_col='x_m',y_col='y_m',n_offsets=n_offsets,
  stratify_methods=TRUE,audit_method_overlap=TRUE),warning=function(w) {
   warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
 })
 observations <- join_soil_main(build_element(d,element),metadata,
  extra_cols=intersect(c(groups,'x_m','y_m','project_name','depth','mesh_size','sample_zone','horizon'),names(metadata)))
 stopifnot(nrow(observations)==nrow(d),sum(result$estimates$n_input)==nrow(d),
  sum(result$method_overlap_support$n_input)==nrow(d),
  sum(result$methods$n_input)==nrow(d),all(is.na(result$declustering$ci_low)),
  all(is.na(result$declustering$ci_high)),all(is.na(result$method_overlap$jaccard_mean[
   result$method_overlap$n_origins_evaluated==0])))
 if(nrow(observations))stopifnot(identical(observations$amethod,d$amethod),
  identical(observations$source_record_number,d$source_record_number))
 value <- list(element=element,signature=signature,result=result,observations=observations,
  warnings=unique(warnings),elapsed_seconds=unname(proc.time()[3]-start),
  status=if(nrow(d))'observed_species_processed' else 'absent_from_selected_cohort')
 temporary <- tempfile(pattern=sprintf('%03d_',i),tmpdir=file.path(folder,'checkpoints'),fileext='.tmp')
 saveRDS(value,temporary);stopifnot(!file.exists(checkpoint),file.rename(temporary,checkpoint))
 cat(stage,medium,element,'complete in',round(value$elapsed_seconds,2),'seconds\n');flush.console()
}
checkpoints <- lapply(seq_along(elements),function(i)readRDS(file.path(folder,'checkpoints',sprintf('%03d.rds',i))))
components <- unique(unlist(lapply(checkpoints,function(x)names(x$result))))
for(name in components) {
 table <- dplyr::bind_rows(lapply(checkpoints,function(x)x$result[[name]]))
 if(ncol(table))fwrite(table,file.path(folder,paste0(name,'.csv')),na='NA')
}
observations <- dplyr::bind_rows(lapply(checkpoints,`[[`,'observations'))
fwrite(observations,file.path(folder,'observations.csv'),na='NA')
status <- data.frame(element=elements,n_input=vapply(checkpoints,function(x)nrow(x$observations),integer(1)),
 elapsed_seconds=vapply(checkpoints,`[[`,numeric(1),'elapsed_seconds'),
 status=vapply(checkpoints,`[[`,character(1),'status'),
 warnings=vapply(checkpoints,function(x)paste(x$warnings,collapse=' | '),character(1)))
fwrite(status,file.path(folder,'run_status.csv'),na='NA')
saveRDS(lapply(checkpoints,`[[`,'result'),file.path(folder,'analysis.rds'))
stopifnot(sum(status$n_input)==nrow(observations),
 sum(status$n_input)==sum(chemistry$species%in%elements))
verify_manifest(source_manifest);verify_manifest(input_manifest)
writeLines(c('Completed method-preserving observational workflow.',
 'All selected source records accounted for; weighted routine confidence limits remain NA.',
 'Reporting design/geological population/site identity remain unestablished; no spatial interval validation.'),
 file.path(folder,'COMPLETED.txt'))
cat('FINISHED',stage,medium,'species',length(elements),'records',nrow(observations),'\n')
