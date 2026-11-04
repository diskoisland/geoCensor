#!/usr/bin/env Rscript
# Research inventory only: no chemistry values, fitted models or spatial tests.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1L] else 'output/usgs_method_domains_2026-11-02'
audit_medium <- if(length(args)>1L) args[2L] else 'soil'
stopifnot(length(args)<=2L,audit_medium%in%c('soil','sed'))
protocol <- file.path(out,'PROTOCOL.md')
source_file <- 'scripts/analysis/audit_usgs_method_domains.R'
stopifnot(file.exists(protocol))
inputs <- file.path('..','data',paste0('ngdb',audit_medium),c('main.csv','bestvalue.csv'))
signature <- list(medium=audit_medium,inputs=tools::md5sum(inputs),
 source=tools::md5sum(source_file),protocol=tools::md5sum(protocol),R=R.version.string,
 data_table=as.character(packageVersion('data.table')))
stopifnot(!anyNA(signature$inputs))
folder <- file.path(out,audit_medium)
manifest <- file.path(folder,'manifest.rds')
if(file.exists(manifest)) {
 stopifnot(identical(readRDS(manifest),signature))
 if(file.exists(file.path(folder,'COMPLETED.txt'))) {
  cat(audit_medium,'completed inventory matches source/input/runtime; no results changed.\n')
  quit(status=0L)
 }
} else {
 dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 stopifnot(!length(list.files(folder,all.files=TRUE,no..=TRUE)))
 saveRDS(signature,manifest)
 file.copy(source_file,file.path(folder,basename(source_file)))
 write.csv(data.frame(file=c(inputs,source_file,protocol),md5=unname(c(signature$inputs,signature$source,signature$protocol))),
  file.path(folder,'source_input_manifest.csv'),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),file.path(folder,'sessionInfo.txt'))
}
normalize_label <- function(x) {
 x <- trimws(as.character(x));x[is.na(x)] <- '';x
}
support_fields <- c('primary_class','depth','mesh_size','prep','collection_protocol',
 'method_collected','sample_source','coord_prec_est','date_collect','date_submitted',
 'geologic_age','stratigraphy','mineralization')
if(audit_medium=='soil') support_fields <- c(support_fields,'sample_zone','horizon')
metadata_fields <- c('lab_id','country','state','lat_wgs84','long_wgs84',
 'project_name','previous_lab_id',support_fields)
cat(audit_medium,'reading sample metadata\n');flush.console()
meta <- fread(inputs[1L],select=metadata_fields,na.strings=NULL,
 colClasses=c(lab_id='character',previous_lab_id='character'))
metadata_input_n <- nrow(meta)
label_fields <- setdiff(metadata_fields,c('lat_wgs84','long_wgs84'))
for(field_name in label_fields) set(meta,j=field_name,value=normalize_label(meta[[field_name]]))
stopifnot(all(nzchar(meta$lab_id)),!anyDuplicated(meta$lab_id))
meta[,is_us:=toupper(country)%in%c('US','USA','UNITED STATES','U.S.')]
for(field_name in c('lat_wgs84','long_wgs84'))
 set(meta,j=field_name,value=suppressWarnings(as.numeric(meta[[field_name]])))
meta[,valid_coordinate:=is.finite(lat_wgs84)&abs(lat_wgs84)<=90&
 is.finite(long_wgs84)&abs(long_wgs84)<=180]
meta[,coordinate_group:=NA_integer_]
meta[is_us&valid_coordinate,coordinate_group:=.GRP,by=.(lat_wgs84,long_wgs84)]
us_meta <- meta[is_us==TRUE]
field_coverage <- rbindlist(lapply(c('project_name','previous_lab_id',support_fields),function(field_name) {
 values <- us_meta[[field_name]]
 data.table(medium=audit_medium,field=field_name,n_samples=nrow(us_meta),
  n_blank=sum(!nzchar(values)),n_literal_unknown=sum(tolower(values)=='unknown'),
  n_nonblank=sum(nzchar(values)),n_distinct_nonblank=uniqueN(values[nzchar(values)]))
}))
fwrite(field_coverage,file.path(folder,'metadata_field_coverage.csv'))
combination_fields <- c('project_name','primary_class','depth','mesh_size','method_collected')
if(audit_medium=='soil') combination_fields <- c(combination_fields,'sample_zone')
combinations <- us_meta[,.(n_samples=.N,n_valid_coordinate_records=sum(valid_coordinate),
 n_exact_coordinate_groups=uniqueN(coordinate_group[!is.na(coordinate_group)])),by=combination_fields]
setorderv(combinations,combination_fields)
fwrite(combinations,file.path(folder,'metadata_support_combinations.csv'))
metadata_projects <- us_meta[,.(n_samples=.N,n_valid_coordinate_records=sum(valid_coordinate),
 n_exact_coordinate_groups=uniqueN(coordinate_group[!is.na(coordinate_group)]),
 n_nonblank_states=uniqueN(state[nzchar(state)])),by=project_name]
setorder(metadata_projects,project_name)
fwrite(metadata_projects,file.path(folder,'metadata_projects.csv'))
cat(audit_medium,'reading best-value identifiers/methods only\n');flush.console()
chem <- fread(inputs[2L],select=c('lab_id','species','amethod','unit'),
 colClasses=c(lab_id='character'),na.strings=NULL)
chem_input_n <- nrow(chem)
for(field_name in names(chem)) set(chem,j=field_name,value=normalize_label(chem[[field_name]]))
meta_index <- match(chem$lab_id,meta$lab_id)
missing_ids <- sum(!nzchar(chem$lab_id))
unmatched <- sum(is.na(meta_index))
matched_us <- !is.na(meta_index)&meta$is_us[meta_index]
matched_non_us <- sum(!is.na(meta_index)&!meta$is_us[meta_index])
stopifnot(!anyNA(matched_us),sum(matched_us)+matched_non_us+unmatched==chem_input_n)
chem[,metadata_index:=meta_index]
chem <- chem[matched_us]
retained <- nrow(chem)
setkey(chem,species)
species_list <- sort(unique(chem$species))
rows <- vector('list',length(species_list))
for(element_index in seq_along(species_list)) {
 audit_species <- species_list[element_index]
 d <- copy(chem[.(audit_species)])
 stopifnot(!anyDuplicated(d$lab_id))
 for(field_name in c('project_name','state','previous_lab_id','coordinate_group',support_fields))
  set(d,j=field_name,value=meta[[field_name]][d$metadata_index])
 domains <- d[,{
  loc <- table(coordinate_group[!is.na(coordinate_group)])
  result <- list(n_records=.N,n_sample_ids=uniqueN(lab_id),
   n_valid_coordinate_records=sum(!is.na(coordinate_group)),
   n_exact_coordinate_groups=length(loc),
   n_repeated_coordinate_groups=sum(loc>1L),
   n_records_at_repeated_coordinates=sum(loc[loc>1L]),
   n_nonblank_states=uniqueN(state[nzchar(state)]),n_blank_state=sum(!nzchar(state)),
   n_resubmission_links=sum(nzchar(previous_lab_id)))
  for(profile_field in support_fields) {
   values <- .SD[[profile_field]]
   result[[paste0('n_blank_',profile_field)]] <- sum(!nzchar(values))
   result[[paste0('n_literal_unknown_',profile_field)]] <- sum(tolower(values)=='unknown')
   result[[paste0('n_distinct_nonblank_',profile_field)]] <- uniqueN(values[nzchar(values)])
  }
  result
 },by=.(species,amethod,unit,project_name)]
 rows[[element_index]] <- domains
 if(element_index%%10L==0L||element_index==length(species_list)) {
  cat(audit_medium,element_index,'/',length(species_list),'species inventoried\n');flush.console()
 }
}
domains <- rbindlist(rows)
setorderv(domains,c('species','amethod','unit','project_name'))
domains[,`:=`(medium=audit_medium,project_recorded=nzchar(project_name),
 method_recorded=nzchar(amethod),unit_recorded=nzchar(unit),
 reporting_design_status='unestablished',physical_site_map_status='unestablished',
 exchangeability_status='requires_geology_and_sampling_review',p_value=NA_real_)]
stopifnot(sum(domains$n_records)==retained,all(domains$n_records==domains$n_sample_ids),
 all(is.na(domains$p_value)),all(domains$n_exact_coordinate_groups<=domains$n_valid_coordinate_records))
fwrite(domains,file.path(folder,'method_project_domains.csv'),na='NA')
method_summary <- domains[,.(n_records=sum(n_records),n_project_labels=.N,
 n_records_blank_project=sum(n_records[!project_recorded])),by=.(species,amethod,unit)]
setorderv(method_summary,c('species','amethod','unit'))
fwrite(method_summary,file.path(folder,'methods_by_species.csv'))
species_summary <- domains[,.(n_records=sum(n_records),n_method_unit_project_groups=.N,
 n_recorded_methods=uniqueN(amethod[nzchar(amethod)]),
 n_recorded_units=uniqueN(unit[nzchar(unit)]),
 n_records_blank_project=sum(n_records[!project_recorded]),
 n_records_blank_method=sum(n_records[!method_recorded]),
 n_records_blank_unit=sum(n_records[!unit_recorded]),
 n_records_blank_mesh=sum(n_blank_mesh_size),n_records_blank_depth=sum(n_blank_depth),
 n_records_blank_protocol=sum(n_blank_collection_protocol)),by=species]
fwrite(species_summary,file.path(folder,'species_summary.csv'))
overview <- data.table(medium=audit_medium,n_input_metadata=metadata_input_n,
 n_us_metadata=nrow(us_meta),n_input_chemistry=chem_input_n,
 n_missing_chemistry_ids=missing_ids,n_unmatched_chemistry=unmatched,
 n_matched_non_us_chemistry=matched_non_us,n_us_chemistry=retained,
 n_species=length(species_list),n_method_unit_project_groups=nrow(domains),
 n_groups_blank_project=sum(!domains$project_recorded),
 n_records_blank_project=sum(domains$n_records[!domains$project_recorded]),
 n_eligible_domains=NA_integer_,p_value=NA_real_)
fwrite(overview,file.path(folder,'overview.csv'),na='NA')
stopifnot(identical(tools::md5sum(inputs),signature$inputs),
 identical(tools::md5sum(source_file),signature$source),
 identical(tools::md5sum(protocol),signature$protocol))
writeLines(c('Method/project inventory complete; all retained records reconciled.',
 'No concentration columns read, no fitted model or p-value computed.',
 'Matching labels and raw nonblank coverage do not establish eligibility.',
 'Reporting designs, physical sites and exchangeability remain unestablished.'),file.path(folder,'COMPLETED.txt'))
print(overview)
