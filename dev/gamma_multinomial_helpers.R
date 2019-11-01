log_ddirichlet <- function(x, alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
}

log_ddirichlet_multinom <- function(y, alpha) {
  alpha_plus <-  sum(alpha);

  lgamma(alpha_plus) + sum(lgamma(alpha + y)) + lgamma(sum(y) + 1) - lgamma(alpha_plus + sum(y)) - sum(lgamma(alpha)) - sum(lgamma(y + 1));
}


plot_integral_approx <- function(N, phi, lambda, res) {
  # if(length(phi) != 1) {
  #   stop("phi")
  # }

  imp_sampling_log <- array(-1, N)
  simple_log <- array(-1, N)
  N_steps = min(100, N)
  step_size = floor(N/N_steps)
  for(n in 1:N) {
    ## The simple way
    gamma_samples <- rgamma(length(lambda), phi, phi/exp(lambda))
    simple_log[n] <- dmultinom(res, prob = gamma_samples, log = TRUE)

    ## Importance sampling
    # proposal_normalized <- MCMCpack::rdirichlet(1, res + 0.5)[1,]
    # proposal_s_alpha = phi * (sum(exp(lambda)) ^ 2) / sum(exp(2 * lambda))
    # proposal_s_beta = sum(exp(2*lambda)) / (phi * sum(exp(lambda)))
    # proposal_s <- rgamma(1, shape = proposal_s_alpha, rate = proposal_s_beta)
    #
    # log_weight <- sum(dgamma(proposal_normalized * proposal_s, shape = phi, rate = phi/exp(lambda), log = TRUE )) - dgamma(proposal_s, shape = proposal_s_alpha, rate = proposal_s_beta, log = TRUE) -
    #   log_ddirichlet(proposal_normalized, res + 0.5)
    #   #log(MCMCpack::ddirichlet(proposal_normalized, res + 0.5))
    #
    # imp_sampling_log[n] <- log_weight + dmultinom(res, prob = proposal_normalized, log = TRUE)

  }
  integral_simple_log <- array(-1, N_steps)
  integral_imp_sampling_log <- array(-1, N_steps)
  for(step in 1:N_steps) {
    integral_simple_log[step] <- logsumexp(simple_log[1:(step * step_size)]) - log(step * step_size)
    integral_imp_sampling_log[step] <- logsumexp(imp_sampling_log[1:(step * step_size)]) - log(step * step_size)
  }
  cat("Final value: ", integral_simple_log[N_steps])
  tibble(n = (1:N_steps) * step_size, integral_simple_log) %>%
    gather("type","value", -n) %>%
    ggplot(aes(x=n, y = value, color = type)) + geom_line()
}

plot_integral_approx_d <- function(N, d) {
  lambda = rnorm(d, 0, 2)
  plot_integral_approx(N,
                       phi = 1 / sqrt(abs(rnorm(1, 0, 1))),
                       lambda = lambda,
                       res = rmultinom(1, size = 1 * d, prob = exp(lambda))[,1]
  )
}
