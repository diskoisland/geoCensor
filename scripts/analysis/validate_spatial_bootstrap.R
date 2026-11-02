#!/usr/bin/env Rscript
# Pilot operating-characteristic study; not a validation of all NGDB designs.
suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE))
folder <- 'output/spatial_bootstrap_validation_2026-10-31'
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
n_sim <- 100L
n_boot <- 199L
settings <- data.frame(
  scenario = c('independent_regular', 'independent_clustered', 'preferential_clustered',
               'correlated_small_blocks', 'correlated_larger_blocks'),
  correlated = c(FALSE, FALSE, FALSE, TRUE, TRUE),
  design = c('regular', 'clustered', 'preferential', 'clustered', 'clustered'),
  block_size = c(1, 1, 1, 1, 2)
)
readr::write_csv(settings, file.path(folder, 'settings.csv'))
writeLines(capture.output(sessionInfo()), file.path(folder, 'session_info.txt'))
writeLines(c(
  'Pilot study: 100 independently generated datasets per setting; 199 bootstrap replicates each.',
  'Target: superpopulation marginal lognormal geometric mean 50 (not a realized finite-domain mean).',
  'Spatial support: 8x8 cells, coordinate units arbitrary; all 64 sites observed.',
  'Within-site correlation: shared Gaussian log-concentration effect sd 0.55 plus independent noise sd 0.2.',
  'Regular design: four observations per site; clustered: 1,4,12 observations sampled independently of grade.',
  'Preferential design: 12 observations where latent effect >0, otherwise two; all sites still observed.',
  'Correlated settings: exponential covariance exp(-distance/0.8) for the site effects.',
  'Cell width 1; bootstrap widths 1 and 2; fixed origin (0,0); two offsets per axis for this pilot.',
  'Repeated site coordinates make declustering invariant to grid offsets in these scenarios.',
  'Left limit exp(log(50)-0.3), upper limit exp(log(50)+1); both censoring directions retained.',
  'The two correlated settings use identical simulated datasets and bootstrap seeds for paired comparison.',
  'Coverage Monte Carlo intervals describe simulation uncertainty, not confidence in a particular real estimate.',
  'Do not infer universal validity or calibrate/tune a method to this study alone.'
), file.path(folder, 'design_notes.txt'))
grid <- expand.grid(x = 0:7, y = 0:7)
distance <- as.matrix(stats::dist(grid))
chol_cov <- chol(0.55^2 * exp(-distance / 0.8) + diag(1e-10, nrow(grid)))
rows <- list()
for (sim in seq_len(n_sim)) {
  checkpoint <- file.path(folder, sprintf('replicate_%03d.rds', sim))
  if (file.exists(checkpoint)) {rows[[sim]] <- readRDS(checkpoint); next}
  result <- list()
  for (j in seq_len(nrow(settings))) {
    setting <- settings[j, ]
    # Paired correlated settings differ only in block size.
    seed_id <- if (j == 5) 4 else j
    set.seed(120000 + sim * 10 + seed_id)
    latent <- if (setting$correlated) as.vector(t(chol_cov) %*% rnorm(64)) else rnorm(64, sd = 0.55)
    counts <- switch(setting$design, regular = rep(4,64),
                      clustered = sample(c(1,4,12),64,replace=TRUE),
                      preferential = ifelse(latent > 0,12,2))
    id <- rep(seq_len(64), counts)
    value <- exp(log(50) + latent[id] + rnorm(length(id),sd=0.2))
    lod <- exp(log(50)-0.3); upper <- exp(log(50)+1)
    x <- data.frame(x = grid$x[id]+0.2, y = grid$y[id]+0.2,
                     detected_val = ifelse(value < lod | value > upper,NA_real_,value),
                     left_cen = value < lod, left_limit = ifelse(value < lod,lod,NA_real_),
                     right_cen = value > upper, right_limit = ifelse(value > upper,upper,NA_real_))
    out <- spatial_block_bootstrap(x,'x','y',cell_size=1,block_size=setting$block_size,
                                    n_boot=n_boot,n_offsets=2,seed=900000+sim*10+seed_id)
    s <- out$summary
    valid <- s$status == 'experimental_interval'
    result[[j]] <- data.frame(simulation=sim,scenario=setting$scenario,block_size=setting$block_size,
                               n_total=nrow(x),fraction_left=mean(x$left_cen),fraction_right=mean(x$right_cen),
                               geomean=s$geomean,ci_low=s$ci_low,ci_high=s$ci_high,n_blocks=s$n_blocks,
                               n_success=s$n_success,status=s$status,
                               covered=valid && s$ci_low <= 50 && s$ci_high >= 50)
  }
  rows[[sim]] <- dplyr::bind_rows(result)
  saveRDS(rows[[sim]],checkpoint)
  if (sim %% 10 == 0) {cat('Completed',sim,'of',n_sim,'simulation datasets in all five settings\n');flush.console()}
}
raw <- dplyr::bind_rows(rows)
readr::write_csv(raw,file.path(folder,'simulation_results.csv'))
summary <- dplyr::bind_rows(lapply(split(raw,raw$scenario),function(x) {
  valid <- x$status == 'experimental_interval'
  n_valid <- sum(valid); hit <- sum(x$covered[valid])
  mc <- if(n_valid) stats::binom.test(hit,n_valid)$conf.int else c(NA_real_,NA_real_)
  data.frame(scenario=x$scenario[1],n_datasets=nrow(x),n_intervals=n_valid,
               coverage=if(n_valid) hit/n_valid else NA_real_,
               coverage_mc_low=mc[1],coverage_mc_high=mc[2],
               mean_width=mean(x$ci_high[valid]-x$ci_low[valid]),
               mean_geomean=mean(x$geomean),mean_left_fraction=mean(x$fraction_left),
               mean_right_fraction=mean(x$fraction_right),
               bootstrap_fit_success_fraction=sum(x$n_success)/(nrow(x)*n_boot))
}))
readr::write_csv(summary,file.path(folder,'coverage_summary.csv'))
print(summary)
writeLines('Pilot simulation complete; see design_notes.txt for scope and limitations.',file.path(folder,'COMPLETED.txt'))
