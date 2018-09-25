#' Simulate resampling count data - given counts of various items, each item is given
#' equal probability to be sampled.
#' @param draw A vector of item counts
#' @param k Number of items to be selected
resample_draw = function(draw, k) {
  if(k <= sum(draw)) {
    stop("The number reads to sample (k) has to be larger than the sum of draw")
  }
  G = length(draw)
  positions = c(0, cumsum(draw))
  read_set = array(-1, sum(draw))
  for(g in 1:G) {
    read_set[(positions[g] + 1):positions[g + 1]] = g
  }

  sampling_result = sample(read_set, k)

  counts = array(-1, G)
  for(g in 1:G) {
    counts[g] = sum(sampling_result == g)
  }

  counts
}

#' @param generator_type either "dm" or "nb_resample" (should be equivalent in practice, but this lets us test)
#' @param model_type either "dm" or "nb" which tells the Stan model which method to use (once again should be equivalent but this lets us test numerical and other issues)
generate_data_dirichlet_multinomial <- function(G, sums, lambda_raw_prior,  theta_given = TRUE, alpha_given = TRUE,
                                  alpha_prior_log_mean = log(5000), alpha_prior_log_sd = 1,
                                  theta_prior_log_mean = log(5), theta_prior_log_sd = 1, generator_type = "dm",
                                  model_type = "nb") {
  N = length(sums)
  counts = array(-1, c(N, G))

  lambda_raw = MCMCpack::rdirichlet(1, rep(lambda_raw_prior,G))

  theta = rlnorm(1,theta_prior_log_mean, theta_prior_log_sd)
  alpha = rlnorm(1, alpha_prior_log_mean, alpha_prior_log_sd)

  if(generator_type == "dm") {
    single_draw_func <- single_draw_dirichlet_multinom
  } else if(generator_type == "nb_resample") {
    single_draw_func <- single_draw_nb_fano_resample
  }

  for(n in 1:N) {
    counts[n,] = single_draw_func(lambda_raw * alpha, theta, sums[n])
  }

  data = list(
    observed = list(
      N = N,
      G = G,
      counts = counts,
      theta_data = if(theta_given) { array(theta, 1) } else {numeric(0)},
      theta_given = theta_given,
      alpha_data = if(alpha_given) { array(alpha, 1) } else {numeric(0)},
      alpha_given = alpha_given,
      lambda_raw_prior = lambda_raw_prior,
      variant = if(model_type == "dm") { 0 } else if (model_type == "nb") { 1} else stop("Invalid model_type")
    ),
    true = list(lambda_raw = lambda_raw,
                theta_param = if(!theta_given) { array(theta, 1) } else {numeric(0)},
                alpha_param = if(!alpha_given) { array(alpha, 1) } else {numeric(0)}
    )
  )

  if(theta_given) {
    data$observed$theta_prior_log_mean = numeric(0)
    data$observed$theta_prior_log_sd = numeric(0)
  } else {
    data$observed$theta_prior_log_mean = array(theta_prior_log_mean, 1)
    data$observed$theta_prior_log_sd = array(theta_prior_log_sd, 1)
  }

  if(alpha_given) {
    data$observed$alpha_prior_log_mean = numeric(0)
    data$observed$alpha_prior_log_sd = numeric(0)
  } else {
    data$observed$alpha_prior_log_mean = array(alpha_prior_log_mean, 1)
    data$observed$alpha_prior_log_sd = array(alpha_prior_log_sd, 1)
  }

  data

}


#'
single_draw_dirichlet_multinom <- function(lambda, k, theta) {
  dirichlet_draw = MCMCpack::rdirichlet(1,  (1/theta) * lambda)
  rmultinom(1, k, dirichlet_draw)
}


single_draw_nb_fano_resample <- function(lambda, k, theta) {
  G = length(lambda)
  for(i in 1:100) {
    raw_counts_draw = rnbinom(G, mu = lambda, size = lambda / theta)
    if(sum(raw_counts_draw) >= num_trials) {
      break;
    }
    if(i == 5) {
      warning("Neded over 5 draws")
    }
  }
  if(sum(raw_counts_draw) < num_trials) {
    stop("Could not draw large enough")
  }

  resample_draw(G, raw_counts_draw, num_trials)
}
