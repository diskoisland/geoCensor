#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1L]else"output/spatial_order_signal_2026-11-01"
stopifnot(file.exists(file.path(out,"COMPLETED.txt")))
manifest<-readRDS(file.path(out,"study_manifest.rds"))
stopifnot(identical(tools::md5sum(names(manifest$sources)),manifest$sources),
  identical(tools::md5sum(names(manifest$inputs)),manifest$inputs))
design<-read.csv(file.path(out,"design_manifest.csv"),stringsAsFactors=FALSE)
stopifnot(identical(unname(tools::md5sum(design$file)),design$md5))
x<-rbind(read.csv(file.path(out,"development_results.csv"),stringsAsFactors=FALSE),
  read.csv(file.path(out,"confirmation_results.csv"),stringsAsFactors=FALSE))
stopifnot(nrow(x)==4000L,length(unique(x$layout_key))==400L,!anyDuplicated(x$key),
  all(is.na(x$p_value[x$status!="available"])),all(x$p_value[x$status=="available"]>=.001))
write.csv(x,file.path(out,"all_results.csv"),row.names=FALSE)
groups<-split(x,interaction(x$phase,x$regime,x$scenario,x$censoring,drop=TRUE))
summary<-do.call(rbind,lapply(groups,function(z) {
  stopifnot(!anyDuplicated(z$layout_key))
  n<-nrow(z);available<-z$status=="available";r<-sum(z$reject%in%TRUE)
  eligible<-all(z$eligible_design)
  ci<-if(eligible)binom.test(r,n)$conf.int else c(NA_real_,NA_real_)
  conditional_ci<-if(eligible&&any(available))binom.test(r,sum(available))$conf.int else c(NA_real_,NA_real_)
  data.frame(phase=z$phase[1L],regime=z$regime[1L],scenario=z$scenario[1L],censoring=z$censoring[1L],
    eligible_design=eligible,null_scenario=z$null_scenario[1L],attempted=n,available=sum(available),
    unavailable=sum(!available),rejections=if(eligible)r else NA_integer_,
    unconditional_rejection=if(eligible)r/n else NA_real_,mc_lower=ci[1L],mc_upper=ci[2L],
    conditional_rejection=if(eligible&&any(available))r/sum(available) else NA_real_,
    conditional_mc_lower=conditional_ci[1L],conditional_mc_upper=conditional_ci[2L],
    no_order_score_variation=sum(z$status=="no_order_score_variation"),
    reporting_guard_refusals=sum(z$status=="heterogeneous_reporting_design"))
}))
rownames(summary)<-NULL
write.csv(summary,file.path(out,"rejection_and_availability.csv"),row.names=FALSE)
counts<-as.data.frame(with(x,table(phase,scenario,censoring,status)),stringsAsFactors=FALSE)
counts<-counts[counts$Freq>0,];names(counts)[5L]<-"attempts"
write.csv(counts,file.path(out,"status_counts.csv"),row.names=FALSE)
null_rows<-x[x$eligible_design&x$null_scenario,]
null_rows$operational_reject<-as.integer(null_rows$reject%in%TRUE)
layout_null<-aggregate(operational_reject~phase+regime+layout_key,data=null_rows,FUN=mean)
names(layout_null)[4L]<-"null_rejection_fraction"
write.csv(layout_null,file.path(out,"layout_null_fractions.csv"),row.names=FALSE)
pooled<-do.call(rbind,lapply(split(layout_null,interaction(layout_null$phase,layout_null$regime)),function(z)
 data.frame(phase=z$phase[1L],regime=z$regime[1L],independent_layouts=nrow(z),
  null_tests=4L*nrow(z),rejection_rate=mean(z$null_rejection_fraction),
  layout_mc_se=sd(z$null_rejection_fraction)/sqrt(nrow(z)))))
write.csv(pooled,file.path(out,"pooled_null_layout_summary.csv"),row.names=FALSE)

confirmation<-summary[summary$phase=="confirmation"&summary$eligible_design,]
confirmation$scenario<-factor(confirmation$scenario,levels=c("independent","site_only","spatial_short","spatial_long"))
confirmation$regime<-factor(confirmation$regime,levels=c("uniform","clustered"))
confirmation$censoring<-factor(confirmation$censoring,levels=c("moderate","heavy"))
plot_rates<-function(data,title) {
 ggplot2::ggplot(data,ggplot2::aes(.data$scenario,.data$unconditional_rejection,color=.data$censoring))+
  ggplot2::geom_hline(yintercept=.05,linetype="dashed",color="grey45")+
  ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$mc_lower,ymax=.data$mc_upper),
    width=.12,position=ggplot2::position_dodge(width=.35))+
  ggplot2::geom_point(size=2.5,position=ggplot2::position_dodge(width=.35))+
  ggplot2::facet_wrap(ggplot2::vars(regime),ncol=2)+
  ggplot2::scale_color_manual(values=c(moderate="#2166ac",heavy="#b35806"),name="Censoring")+
  ggplot2::scale_y_continuous(labels=function(z)paste0(format(100*z,trim=TRUE,scientific=FALSE),"%"))+
  ggplot2::labs(x="Generating scenario",y="Operational rejection rate",title=title,
    subtitle="Confirmation only: 150 independent layouts per regime/scenario/severity; nominal level 5%.",
    caption="Bars are 95% binomial Monte Carlo intervals for simulation rates, not geochemical confidence limits.\nUnavailable tests remain NA and count as no operational rejection; availability is reported separately.")+
  ggplot2::theme_minimal(base_size=11)+ggplot2::theme(legend.position="bottom")
}
null_plot<-plot_rates(confirmation[confirmation$null_scenario,],
  "Fresh confirmation: false rejection under exchangeable site-block nulls")
power_plot<-plot_rates(confirmation[!confirmation$null_scenario,],
  "Fresh confirmation: detecting positive spatial ordering")
ggplot2::ggsave(file.path(out,"confirmation_null_rates.png"),null_plot,width=11,height=6.5,dpi=180)
ggplot2::ggsave(file.path(out,"confirmation_power.png"),power_plot,width=11,height=6.5,dpi=180)
saveRDS(list(null=null_plot,power=power_plot),file.path(out,"plots.rds"))
fmt<-function(z)if(is.finite(z))sprintf("%.1f%%",100*z)else"NA"
table<-vapply(seq_len(nrow(confirmation)),function(i) {
 z<-confirmation[i,]
 sprintf("| %s | %s | %s | %d/%d | %s | %s–%s |",z$regime,z$scenario,z$censoring,z$available,
   z$attempted,fmt(z$unconditional_rejection),fmt(z$mc_lower),fmt(z$mc_upper))
},"")
writeLines(c("# Censor-aware site-block spatial-order evidence — November 1, 2026","",
 "This research test asks whether observed site orderings show positive spatial concordance under a homogeneous, exchangeable whole-site null. It fits no concentration model or covariance range and substitutes no censored concentration. An exact value is scored above another interval only when its bounds establish that strict ordering; overlap/ties contribute no definite ordering. Whole paired site vectors are preserved.","",
 "## Eligibility and interpretation","",
 "The null requires exchangeable site observation vectors conditional on fixed exogenous coordinates and reporting design. Correlated repeats are allowed within sites. Equal covariances alone, matching metadata or a software assertion do not establish exchangeability. Geology/method mixtures, trends, response-selected locations or reporting mechanisms can invalidate the interpretation. Every record's reporting design, including exact records, must be supplied. Heterogeneous reporting designs and missing inputs return NA; sites are not silently dropped.","",
 "Read the [scientific interpretation and limits](INTERPRETATION.md) alongside these numerical results.","",
 "A rejection concerns this positive spatial-ordering statistic under the stated null. It does not prove an exponential covariance family, determine a range, distinguish trend from residual covariance, correct preferential sampling, validate a geometric mean or supply its confidence limits. Non-rejection does not establish independence; heavy censoring/ties can leave little observable contrast. All-left/constant-order scores produce an unavailable test, not p=1. Censoring bounds can be informative even with no exact records, independently of the public concentration-fit count rule.","",
 "## Fresh development and confirmation","",
 "All 400 layouts/permutation plans and 2,000 response fields were generated and fingerprinted before outcomes. Development used 50 layouts per geometry regime; confirmation used 150 per regime, with no tuning or source change. Each layout has 24 sites/two records per site and five scenarios, paired across censoring severity. Four thousand attempts comprise 3,200 eligible-design analyses and 800 deliberate reporting-design refusals, not 4,000 independent surveys. Scenarios/severity partners share layouts and must not be pooled as independent replicates.","",
 sprintf("Observed reporting-design refusals: %d/800. Eligible-design available tests: %d/3200; absent-order-score outcomes: %d. All unavailable p-values remain NA. Confirmation has %d eligible-design analyses on 300 independent layouts, with one response per layout in each displayed stratum.",
   sum(x$status=="heterogeneous_reporting_design"),sum(x$eligible_design&x$status=="available"),
   sum(x$status=="no_order_score_variation"),sum(x$phase=="confirmation"&x$eligible_design)),"",
 "| Geometry | Scenario | Censoring | Available | Operational rejection | 95% Monte Carlo interval |",
 "| --- | --- | --- | ---: | ---: | --- |",table,"",
 "The table reports false rejection for independent/site-only nulls and power for positive-range fields. Unavailable p-values are not changed to one: operational summaries count no rejection and retain separate availability/conditional-rate columns. Intervals quantify simulation-rate uncertainty for independent-layout strata, not geochemical parameter uncertainty. There is no across-elements/domains multiplicity claim.","",
 "pooled_null_layout_summary.csv pools the four null tests within each layout before calculating sampling Monte Carlo SE. Those four dependent tests are not counted as independent surveys.","",
 "![Confirmation null rejection](confirmation_null_rates.png)","",
 "![Confirmation power](confirmation_power.png)","",
 "## Calibration and retained evidence","",
 "Three fixed geometry distance scales are combined by the maximum Moran-type score quadratic. Every permutation recomputes that same maximum, calibrating this internal scale search rather than using an unadjusted smallest per-scale p-value. Coordinates stay fixed; complete site vectors move together. With 999 uniform random permutations and the observed identity, p=(1+random exceedances)/(999+1), using conservative ties. P-values cannot be zero. These are restricted null randomizations, not a spatial bootstrap for confidence intervals.","",
 "The group-invariance and identity-inclusion justification follows [Hemerik & Goeman (2018)](https://doi.org/10.1007/s11749-017-0571-1). The nonzero random-permutation calculation is also discussed by [Phipson & Smyth (2010), author preprint](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf). Those results support the machinery under assumptions, not real-NGDB eligibility or power of this new score.","",
 "Thirteen exact-rank/equivariance/full-group/plus-one comparisons and missing/design/RNG guards passed before validation. Exhaustive six-site enumeration checked all 720 assignments at four levels, including ties and the maximum over scales. Every checkpoint retains observed bounds, full reporting design, score vectors and permutation statistics. Frozen source/runtime and design manifests guard workers and resume. Sources and all outcomes/failures are preserved.","",
 "No real USGS test or production API was added. Confidence limits remain experimental. A future application needs substantive domain/exchangeability review and a declared multiplicity strategy across elements/domains. Any method revision informed by these results requires a new prespecified confirmation study."),file.path(out,"README.md"))
invisible(file.copy("scripts/analysis/report_spatial_block_order_test.R",file.path(out,"source_snapshot"),overwrite=TRUE))
cat("Development/confirmation rates, availability, guard refusals and Monte Carlo intervals reported.\n")
