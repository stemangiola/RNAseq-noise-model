sbc <- function(model, generator, N_steps, ...) {
  1:N_steps %>% map_dfr(function(x) {
    data <- generator()
    fit <- rstan::sampling(model, data = data$observed, ...)
    samples <- rstan::extract(fit, pars = names(data$true))
    eval <- evaluate_all_params(samples, data$true)
    eval %>% mutate(run = x)
  })
}

plot_sbc <- function(results, n_bins = 10, caption = NULL) {
  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  n_simulations <- length(unique(results$run))
  CI = qbinom(c(0.005,0.5,0.995), size=n_simulations,prob  =  1/n_bins)
  lower = CI[1]
  mean = CI[2]
  upper = CI[3]

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  print(results %>%
          ggplot(aes(x = order_within / 100)) +
          geom_segment(aes(x=0,y=mean,xend=1,yend=mean),colour="grey25") +
          geom_polygon(data=data.frame(x=c(-0.1,0,-0.1,1.1,1.00,1.1,-0.1),y=c(lower,mean,upper,upper,mean,lower,lower)),aes(x=x,y=y),fill="grey45",color="grey25",alpha=0.5) +
          geom_histogram(breaks=seq(0,1,length.out = n_bins + 1) ,fill="#A25050",colour="black") +
          facet_wrap(~param_name, scales = "free_y") +
          ggtitle("Posterior order within 100 samples")
  )

  point_alpha <- 1 / ((n_simulations * 0.03) + 1)
  print(results %>%
          ggplot(aes(x = true_value, y = median)) + geom_point(alpha = point_alpha) +
          geom_abline(slope = 1, intercept = 0, color = "blue") +
          facet_wrap(~param_name, scales = "free")  +
          ggtitle(paste0("Median of marginal posteriors vs. true value - ", caption))
  )
}
