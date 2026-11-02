#!/usr/bin/env Rscript
source("scripts/analysis/spatial_block_order_test_reference.R")
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "output/spatial_order_signal_2026-11-01"
mode <- if (length(args)>=2L) args[2L] else "init"
phase <- if (length(args)>=3L) args[3L] else "development"
if (length(args)>3L || !mode %in% c("init","worker1","worker2","finish") ||
    !phase %in% c("development","confirmation")) stop("Invalid study mode/phase")
RNGkind("Mersenne-Twister","Inversion","Rejection")
config <- list(sites=24L,repeats=2L,regimes=c("uniform","clustered"),
  phases=c(development=50L,confirmation=150L),workers=2L,
  scenarios=c("independent","site_only","spatial_short","spatial_long","varying_reporting"),
  truth=log(50),sigma=sqrt(.55^2+.2^2),nugget=.2^2/(.55^2+.2^2),
  offsets=c(moderate=-.3,heavy=.6),upper_offset=1,scales=c(.05,.15,.4),
  permutations=999L,alpha=.05)
sources <- c("scripts/analysis/validate_spatial_block_order_test.R",
  "scripts/analysis/spatial_block_order_test_reference.R",
  "scripts/analysis/verify_spatial_block_order_test_reference.R")
inputs <- file.path(out,c("PROTOCOL.md","verification.csv","VERIFIED.txt","exhaustive_group_audit.csv"))
stopifnot(all(file.exists(inputs)))
signature <- list(config=config,sources=tools::md5sum(sources),inputs=tools::md5sum(inputs),
  R=R.version.string,platform=R.version$platform,RNG=RNGkind())
manifest <- file.path(out,"study_manifest.rds")
if(file.exists(manifest)) {
  if(!identical(readRDS(manifest),signature))stop("Order-test source/protocol/runtime changed")
} else {
  if(mode!="init" || length(setdiff(list.files(out),basename(inputs))))stop("Unverified study output")
  dir.create(file.path(out,"source_snapshot"))
  invisible(file.copy(sources,file.path(out,"source_snapshot")))
  saveRDS(signature,manifest)
  write.csv(data.frame(file=c(sources,inputs),md5=unname(c(signature$sources,signature$inputs))),
    file.path(out,"source_input_manifest.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
}
verify_design <- function() {
  inventory<-read.csv(file.path(out,"design_manifest.csv"),stringsAsFactors=FALSE)
  stopifnot(identical(unname(tools::md5sum(inventory$file)),inventory$md5))
}
if(file.exists(file.path(out,"INITIALIZED.txt")))verify_design()
if(file.exists(file.path(out,"COMPLETED.txt"))) {
  cat("Completed order-test study matches source/runtime/design; no results changed.\n");quit(status=0L)
}
for(d in c("design","checkpoints"))dir.create(file.path(out,d),showWarnings=FALSE)
if(mode=="init") {
  if(file.exists(file.path(out,"INITIALIZED.txt")))stop("Already initialized")
  jobs<-list();j<-0L
  for(p in seq_along(config$phases))for(r in seq_along(config$regimes))for(i in seq_len(config$phases[p])) {
    name<-sprintf("%s_%s_%03d",names(config$phases)[p],config$regimes[r],i)
    path<-file.path(out,"design",paste0(name,".rds"))
    if(!file.exists(path)) {
      set.seed(910010000L+p*1000000L+r*10000L+i)
      sites<-if(r==1L)matrix(runif(48L),ncol=2L)else
        rbind(matrix(runif(36L,0,.25),ncol=2L),matrix(runif(12L),ncol=2L))
      plan<-block_order_plan_reference(sites,config$scales,config$permutations,
        930010000L+p*1000000L+r*10000L+i)
      id<-rep(1:config$sites,each=config$repeats);xy<-sites[id,];distance<-as.matrix(dist(xy))
      fields<-list()
      for(s in seq_along(config$scenarios)) {
        scenario<-config$scenarios[s]
        kernel<-switch(scenario,independent=diag(length(id)),site_only=outer(id,id,"==")*1,
          spatial_short=exp(-distance/.12),spatial_long=exp(-distance/.3),
          varying_reporting=outer(id,id,"==")*1)
        nu<-if(scenario=="independent")1 else config$nugget
        covariance<-config$sigma^2*((1-nu)*kernel+diag(nu,length(id)))
        seed<-920010000L+p*1000000L+r*10000L+10L*i+s
        set.seed(seed)
        y<-config$truth+as.vector(t(chol(covariance))%*%rnorm(length(id)))
        fields[[scenario]]<-list(log_values=matrix(y,config$sites,config$repeats,byrow=TRUE),
          covariance=covariance,seed=seed)
      }
      saveRDS(list(key=name,phase=names(config$phases)[p],regime=config$regimes[r],
        sites=sites,plan=plan,fields=fields,layout_seed=910010000L+p*1000000L+r*10000L+i),path)
    }
    j<-j+1L;jobs[[j]]<-data.frame(key=name,phase=names(config$phases)[p],regime=config$regimes[r],
      index=i,worker=(i-1L)%%2L+1L)
    if(j%%50L==0L)cat("Initialized",j,"layouts\n")
  }
  write.csv(do.call(rbind,jobs),file.path(out,"jobs.csv"),row.names=FALSE)
  files<-c(file.path(out,"jobs.csv"),list.files(file.path(out,"design"),full.names=TRUE))
  write.csv(data.frame(file=files,md5=unname(tools::md5sum(files))),file.path(out,"design_manifest.csv"),row.names=FALSE)
  writeLines("All400 layouts/plans and2000 response fields frozen before outcomes.",file.path(out,"INITIALIZED.txt"))
  cat("All development/confirmation responses initialized without observing test outcomes.\n");quit(status=0L)
}
stopifnot(file.exists(file.path(out,"INITIALIZED.txt")))
jobs<-read.csv(file.path(out,"jobs.csv"),stringsAsFactors=FALSE)
selected<-jobs[jobs$phase==phase,]
if(phase=="confirmation" && !file.exists(file.path(out,"DEVELOPMENT_COMPLETED.txt")))
  stop("Finish frozen development phase before confirmation")
if(mode=="finish") {
  paths<-file.path(out,"checkpoints",paste0(selected$key,".rds"))
  stopifnot(all(file.exists(paths)))
  checkpoints<-lapply(paths,readRDS);rows<-do.call(rbind,lapply(checkpoints,`[[`,"rows"))
  expected<-nrow(selected)*length(config$scenarios)*length(config$offsets)
  stopifnot(nrow(rows)==expected,!anyDuplicated(rows$key),
    all(is.na(rows$p_value[rows$status!="available"])),
    all(rows$p_value[rows$status=="available"]>0),
    all(rows$status[rows$scenario=="varying_reporting"]=="heterogeneous_reporting_design"))
  write.csv(rows,file.path(out,paste0(phase,"_results.csv")),row.names=FALSE)
  stopifnot(identical(tools::md5sum(sources),signature$sources),identical(tools::md5sum(inputs),signature$inputs))
  verify_design()
  marker<-if(phase=="development")"DEVELOPMENT_COMPLETED.txt" else "COMPLETED.txt"
  writeLines(paste(phase,":",nrow(selected),"layouts",nrow(rows),"attempts; no tuning or data fill."),
    file.path(out,marker))
  cat(phase,"finalized:",nrow(rows),"attempts.\n");quit(status=0L)
}
worker<-as.integer(sub("worker","",mode))
for(j in which(selected$worker==worker)) {
  job<-selected[j,];path<-file.path(out,"checkpoints",paste0(job$key,".rds"))
  if(file.exists(path)){stopifnot(readRDS(path)$key==job$key);next}
  design<-readRDS(file.path(out,"design",paste0(job$key,".rds")))
  results<-list();rows<-list();k<-0L
  for(scenario in config$scenarios)for(severity in names(config$offsets)) {
    y<-design$fields[[scenario]]$log_values
    L<-matrix(config$truth+config$offsets[[severity]],config$sites,config$repeats)
    U<-matrix(config$truth+config$upper_offset,config$sites,config$repeats)
    if(scenario=="varying_reporting")L<-L+.5*(design$sites[,1L]-.5)
    lower<-ifelse(y<L,-Inf,ifelse(y>U,U,y));upper<-ifelse(y>U,Inf,ifelse(y<L,L,y))
    fit<-tryCatch(spatial_block_order_test_reference(lower,upper,design$sites,L,U,
      exchangeable_under_null=TRUE,plan=design$plan),error=function(e)
        list(status="numerical_failure",p_value=NA_real_,statistic=NA_real_,
          definite_cross_site_record_pairs=NA_real_,message=conditionMessage(e)))
    k<-k+1L;key<-paste(job$key,scenario,severity,sep="_")
    results[[k]]<-list(key=key,lower=lower,upper=upper,lower_reporting=L,upper_reporting=U,result=fit)
    available<-fit$status=="available"
    rows[[k]]<-data.frame(key=key,layout_key=job$key,phase=phase,regime=job$regime,
      scenario=scenario,censoring=severity,eligible_design=scenario!="varying_reporting",
      null_scenario=scenario%in%c("independent","site_only"),
      n_sites=config$sites,n_records=length(lower),n_exact=sum(lower==upper),
      n_left=sum(is.infinite(lower)),n_right=sum(is.infinite(upper)),
      status=fit$status,p_value=fit$p_value,statistic=fit$statistic,
      definite_pairs=fit$definite_cross_site_record_pairs,
      reject=if(available)fit$p_value<=config$alpha else NA,
      message=if(is.null(fit$message))"" else fit$message)
  }
  saveRDS(list(key=job$key,rows=do.call(rbind,rows),results=results),path)
  if(job$index%%10L==0L){cat(phase,job$regime,job$index,"worker",worker,"layouts complete\n");flush.console()}
}
verify_design()
stopifnot(identical(tools::md5sum(sources),signature$sources),identical(tools::md5sum(inputs),signature$inputs))
cat(phase,"worker",worker,"complete.\n")
