#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr);library(ggplot2)})
folder <- 'output/spatial_bootstrap_validation_2026-10-31'
stopifnot(file.exists(file.path(folder,'COMPLETED.txt')))
s <- read_csv(file.path(folder,'coverage_summary.csv'),show_col_types=FALSE)
x <- read_csv(file.path(folder,'simulation_results.csv'),show_col_types=FALSE)
stopifnot(nrow(s)==5,nrow(x)==500,all(s$n_datasets==100),all(s$n_intervals<=100))
labels <- c(independent_regular='Independent sites; regular density',
             independent_clustered='Independent sites; uneven density',
             preferential_clustered='Independent sites; grade-dependent density',
             correlated_small_blocks='Correlated sites; smaller blocks',
             correlated_larger_blocks='Correlated sites; larger blocks')
s$label <- factor(labels[s$scenario],levels=rev(unname(labels)))
s$design <- ifelse(grepl('^correlated',s$scenario),'Correlation between sites','Independent sites')
p <- ggplot(s,aes(y=label,x=coverage,color=design))+
  geom_vline(xintercept=0.95,linetype='dashed',color='grey45')+
  geom_segment(aes(x=coverage_mc_low,xend=coverage_mc_high,yend=label),linewidth=0.8)+
  geom_point(size=3)+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.1),labels=function(v)paste0(round(v*100),'%'))+
  scale_color_manual(values=c('Correlation between sites'='#bb5c24','Independent sites'='#176c93'))+
  labs(title='Coverage of nominal 95% spatial bootstrap intervals',
       subtitle='Pilot simulations with clustered, left- and right-censored lognormal data',
       x='Fraction of intervals containing the known geometric mean',y=NULL,color=NULL,
       caption='100 datasets per setting; 199 bootstrap replicates per dataset. Bars: exact 95% Monte Carlo intervals.\nDashed line: nominal coverage. All 64 sites were observed; unsampled-area bias was not tested.')+
  theme_minimal(base_size=11)+theme(legend.position='bottom',panel.grid.major.y=element_blank())
ggsave(file.path(folder,'coverage.png'),p,width=11,height=5.5,dpi=160)
lines <- c('# Spatial block-bootstrap pilot study','',
  'The experimental procedure is implemented and has passed software tests. This study does not establish general validity for NGDB sampling designs.','',
  'Each setting used 100 independently simulated datasets and 199 resamples per dataset. The target was the superpopulation geometric mean of 50. Both left and right censoring were included.','',
  '| Setting | Intervals returned | Coverage | 95% Monte Carlo interval | Mean width |',
  '| --- | ---: | ---: | --- | ---: |')
for (i in seq_len(nrow(s))) {
  z<-s[i,]
  lines<-c(lines,sprintf('| %s | %d/100 | %.0f%% | %.1f–%.1f%% | %.1f |',
                         labels[z$scenario],z$n_intervals,100*z$coverage,
                         100*z$coverage_mc_low,100*z$coverage_mc_high,z$mean_width))
}
lines<-c(lines,'','![Coverage comparison](coverage.png)','',
  'Interpretation: use the observed coverage, not replicate convergence alone, to judge the procedure. Undercoverage means the intervals are too narrow or otherwise miscalibrated in that setting. Larger blocks may help, but fewer occupied blocks and dependence across boundaries remain concerns.','',
  'Keep this implementation experimental. Do not automatically add these intervals to the full USGS analysis. Next validation should assess genuinely irregular locations, unsampled areas, stronger informative site selection, block-origin sensitivity, additional correlation scales, and more Monte Carlo replications. Compare a design-appropriate bootstrap or spatial model rather than tuning solely to this pilot.','',
  'Algorithm: draw occupied nonoverlapping blocks with replacement; recompute declustering weights on unique selected rows; multiply by block draw counts; refit the weighted censored model; take percentile limits on the log scale. Duplicate block draws are not cancelled by declustering repeated coordinates.','',
  'The independent-site settings still have shared effects among repeated observations at the same site. The correlated settings add exponential spatial covariance between sites and compare block widths 1 and 2 using the same datasets. Declustering cell width is 1.','',
  'Preferential density changes the observation count at already observed sites; it does not simulate preferential omission of whole sites. All sites lie on an 8-by-8 grid. These restrictions make this an initial diagnostic study.','',
  'Code: scripts/analysis/validate_spatial_bootstrap.R and scripts/analysis/report_spatial_bootstrap.R. Dataset and bootstrap seeds are fixed in the simulation script. Checkpoints require a new output directory if the design or implementation changes.','',
  'Files: settings.csv, design_notes.txt, simulation_results.csv, coverage_summary.csv, coverage.png, session_info.txt, and dataset checkpoints.')
writeLines(lines,file.path(folder,'README.md'))
cat('Simulation output checks and report completed.\n')
