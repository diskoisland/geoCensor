constant_coded_sample <- function() {
  tibble::tibble(detected_val=rep(10,8),left_cen=FALSE,left_limit=NA_real_,
    right_cen=FALSE,right_limit=NA_real_,amethod='ES_SQ')
}

test_that('constant coding never produces a positive-scale lognormal fit', {
  x <- constant_coded_sample()
  for(w in list(rep(1,8),c(.01,.02,.03,.04,.1,.2,.25,.35),seq_len(8))) {
    x$w <- w
    fit <- fit_cen_lognorm(x,weight_col='w')
    expect_equal(fit$method,'failed_fit')
    expect_true(all(is.na(fit[c('geomean','ci_low','ci_high')])))
    expect_equal(fit$n_total,8)
    expect_equal(fit$n_semiquantitative_used,8)
    expect_true(fit$review_semiquantitative)
  }
  expect_equal(fit_cen_lognorm(x)$method,'failed_fit')
  q <- fit_lognorm_with_quants(x)
  expect_equal(q$method,'failed_fit')
  expect_true(all(is.na(q[c('meanlog','sdlog','geomean','ci_low','ci_high','p10','median','p90','p95','p99')])))
  expect_error(plot_lognormal_fit(x),'failed fit')
})

test_that('mixed censoring is rejected only when all intervals include the point', {
  x <- constant_coded_sample()
  x$detected_val[5:8] <- NA_real_
  x$left_cen[5:6] <- TRUE;x$left_limit[5:6] <- c(10,20)
  x$right_cen[7:8] <- TRUE;x$right_limit[7:8] <- c(10,5)
  fit <- fit_cen_lognorm(x)
  expect_equal(fit$method,'failed_fit')
  expect_equal(c(fit$n_total,fit$n_left,fit$n_right),c(8,2,2))
  x$left_limit[5] <- 5;x$right_limit[7] <- 20
  fit <- fit_cen_lognorm(x)
  expect_equal(fit$method,'survreg_lognormal_MLE_L+R')
  expect_true(is.finite(fit$geomean))
})

test_that('zero weights and distinct supplied values determine the active boundary', {
  x <- dplyr::bind_rows(constant_coded_sample(),constant_coded_sample()[1,])
  x$detected_val[9] <- 100;x$w <- c(rep(1,8),0)
  expect_equal(fit_cen_lognorm(x,weight_col='w')$method,'failed_fit')
  x$w[9] <- 1
  expect_equal(fit_cen_lognorm(x,weight_col='w')$method,'weighted_survreg_lognormal')
  expect_false(geoCensor:::.censored_zero_scale_boundary(c(10,10+1e-8),c(10,10+1e-8)))
})

test_that('boundary AIC fits are excluded without removing fixed-scale families', {
  x <- constant_coded_sample()
  result <- compare_censored_distributions(x,c('lognormal','loggaussian','weibull',
    'loglogistic','gaussian','logistic','extreme','t','exponential','rayleigh'))
  estimated <- !result$distribution%in%c('exponential','rayleigh')
  expect_true(all(!result$converged[estimated]))
  expect_true(all(is.na(result$AIC[estimated])))
  expect_true(all(is.na(result$akaike_weight[estimated])))
  expect_true(all(!result$best[estimated]))
  expect_true(all(result$message[estimated]=='zero_scale_boundary'))
  expect_true(all(result$converged[!estimated]))
  expect_true(all(is.finite(result$AIC[!estimated])))
  expect_equal(result$n_semiquantitative_used,rep(8,nrow(result)))
})

test_that('empirical bounds retain constant records when a model is unavailable', {
  x <- constant_coded_sample()
  empirical <- censored_cdf_bounds(x,c(9,10,11))
  expect_equal(empirical$cdf$cdf_lower,c(0,1,1))
  expect_equal(empirical$cdf$cdf_upper,c(0,1,1))
  expect_equal(nrow(x),8)
})
