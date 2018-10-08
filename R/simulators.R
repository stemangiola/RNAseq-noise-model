

#' @param generator_type either "dm" or "nb_resample" (should be equivalent in practice, but this lets us test)
#' @param model_type either "dm" or "nb" which tells the Stan model which method to use (once again should be equivalent but this lets us test numerical and other issues)
generate_data_dirichlet_multinomial_base <- function(G, sums, lambda_raw_prior,  theta_given = TRUE, alpha_given = TRUE,
                                  alpha_prior_log_mean = 3, alpha_prior_log_sd = 1,
                                  theta_prior_log_mean = log(5), theta_prior_log_sd = 1, generator_type = "dm",
                                  model_type = "nb") {
  N = length(sums)
  counts = array(-1, c(N, G))

  #lambda_raw = MCMCpack::rdirichlet(1, rep(lambda_raw_prior,G))
  lambda_raw = rlnorm(G, 0, lambda_raw_prior)
  lambda_raw = lambda_raw / sum(lambda_raw)

  theta = rlnorm(1,theta_prior_log_mean, theta_prior_log_sd)
  alpha = rlnorm(1, alpha_prior_log_mean, alpha_prior_log_sd)

  if(generator_type == "dm") {
    single_draw_func <- single_draw_dirichlet_multinom
  } else if(generator_type == "nb_resample") {
    single_draw_func <- single_draw_nb_fano_resample
  } else {
    stop("Invalid generator_type")
  }

  for(n in 1:N) {
    counts[n,] = single_draw_func(lambda_raw * sums[n] * alpha, sums[n], theta)
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


single_draw_nb_fano_resample <- function(lambda, k, theta, num_trials = 100) {
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
  if(sum(raw_counts_draw) < k) {
    stop("Could not draw large enough")
  }

  resample_draw(raw_counts_draw, k)
}

single_draw_nb_resample <- function(lambda, k, phi, num_trials = 100) {
  G = length(lambda)
  for(i in 1:100) {
    raw_counts_draw = rnbinom(G, mu = lambda, size = phi)
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

  resample_draw(raw_counts_draw, k)
}

generate_data_nb_multinomial <- function(G, sums, lambda_raw_sigma,  phi_given = TRUE, alpha_given = TRUE,
                                                     alpha_prior_log_mean = 3, alpha_prior_log_sd = 1,
                                                     phi_prior_log_mean = log(5), phi_prior_log_sd = 1) {
  N = length(sums)
  counts = array(-1, c(N, G))

  #lambda_raw = MCMCpack::rdirichlet(1, rep(lambda_raw_prior,G))
  lambda_raw = softmax(rnorm_sum_to_zero(G) * lambda_raw_sigma)

  phi = rlnorm(1,phi_prior_log_mean, phi_prior_log_sd)
  alpha = 1 + rlnorm(1, alpha_prior_log_mean, alpha_prior_log_sd)

  for(n in 1:N) {
    counts[n,] = single_draw_nb_resample(lambda_raw * sums[n] * alpha, sums[n], phi)
  }

  data = list(
    observed = list(
      N = N,
      G = G,
      counts = counts,
      phi_data = if(phi_given) { array(phi, 1) } else {numeric(0)},
      phi_given = phi_given,
      alpha_data = if(alpha_given) { array(alpha, 1) } else {numeric(0)},
      alpha_given = alpha_given,
      lambda_raw_sigma = lambda_raw_sigma
    ),
    true = list(lambda_raw = lambda_raw,
                phi_param = if(!phi_given) { array(phi, 1) } else {numeric(0)},
                alpha_param = if(!alpha_given) { array(alpha, 1) } else {numeric(0)}
    )
  )

  if(phi_given) {
    data$observed$phi_prior_log_mean = numeric(0)
    data$observed$phi_prior_log_sd = numeric(0)
  } else {
    data$observed$phi_prior_log_mean = array(phi_prior_log_mean, 1)
    data$observed$phi_prior_log_sd = array(phi_prior_log_sd, 1)
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

generate_data_multinomial <- function(G, sums, sigma_prior_type = "cauchy", sigma_prior_sigma = 2, use_lambda_mu = FALSE) {

  N = length(sums)

  if(sigma_prior_type == "cauchy") {
    lambda_sigma <- abs(rcauchy(1, 0, sigma_prior_sigma))
    sigma_prior_type_int = 0
  } else if(sigma_prior_type == "normal") {
    lambda_sigma <- abs(rnorm(1, 0, sigma_prior_sigma))
    sigma_prior_type_int = 1
  } else {
    stop("Invlaid sigma_prior_type")
  }

  lambda <- rnorm_sum_to_zero(G) * lambda_sigma

  counts = array(-1, c(N, G))
  for(n in 1:N) {
    counts[n, ] = rmultinom(1, sums[n], softmax(lambda))
  }

  data = list(
    observed = list(
      N = N,
      G = G,
      counts = counts,
      my_prior = c(0,0),
      omit_data = 0,
      exposure = rowSums(counts),
      sigma_prior_type = sigma_prior_type_int,
      sigma_prior_sigma = sigma_prior_sigma
    ),
    true = list(
      lambda_sigma = lambda_sigma,
      lambda = lambda
    )
  )
}
