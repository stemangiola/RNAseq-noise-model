functions {
  vector lognormal_mean_match(vector y, vector phi, real[] x_r, int[] x_i) {
    vector[1] z;
    z[1] = (exp(square(y[1])) - 1) * exp(2 * phi[1] + square(y[1])) - phi[2];
    return(z);
  }
  
}

data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real<lower=0> lambda_prior;
  int<lower=0,upper=1> phi_given;
  
  vector<lower=0>[phi_given ? 1 : 0] phi_data;
  vector<lower=0>[phi_given ? 0 : 1] phi_prior_log_mean;
  vector<lower=0>[phi_given ? 0 : 1] phi_prior_log_sd;

  int<lower=0,upper=1> alpha_given;
  vector<lower=0>[alpha_given ? 1 : 0] alpha_data;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_mean;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_sd;
}

transformed data {
  vector[1] mean_match_guess;
  real empty_real[0];
  int empty_int[0];
  
  mean_match_guess[1] = 1;//log(alpha + square(alpha / N) * N) - alpha;
}

parameters {
  simplex[G] lambda_raw;
  vector<lower=0>[phi_given ? 0 : 1] phi_param_raw;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_param_raw;
}

transformed parameters {
  vector<lower=0>[phi_given ? 0 : 1] phi_param = exp(phi_param_raw * phi_prior_log_sd + phi_prior_log_mean);
  vector<lower=0>[alpha_given ? 0 : 1] alpha_param = exp(alpha_param_raw * alpha_prior_log_sd + alpha_prior_log_mean);
  real phi = phi_given ? phi_data[1] : phi_param[1];
  real alpha = alpha_given ? alpha_data[1] : alpha_param[1];
  vector[G] lambda = lambda_raw * alpha;
  real<lower=0> phi_bar = phi * square(sum(lambda)) / (sum(square(lambda)));
    
  //Note that alpha = sum(lambda)
  real nb_sum_variance = alpha + sum(square(lambda)) / phi;

  real log_bar_helper = 1 + nb_sum_variance/square(alpha);
  real mean_log_bar = log(alpha / sqrt(log_bar_helper));
  real sd_log_bar = sqrt(log(log_bar_helper));

  /* === Alternative way to determine params, seems to perform a tiny bit better, but is more complicated */
  
  // real mean_log_bar = log(sum(lambda));
  // vector[2] solver_params = to_vector({mean_log_bar, nb_sum_variance});
  // real sd_log_bar = algebra_solver(
  //   lognormal_mean_match,
  //   mean_match_guess,
  //   solver_params,
  //   empty_real,
  //   empty_int)
  //   [1];
}

model {
  lambda_raw ~ dirichlet(rep_vector(lambda_prior, G));
  
  if(!alpha_given) {
    alpha_param_raw ~ normal(0,1);
  }
  
  if(!phi_given) {
    phi_param_raw ~ normal(0, 1);
  }
  
  for(n in 1:N) {
    counts[n,] ~ neg_binomial_2(lambda, phi);
    //target += neg_binomial_2_lpmf(counts[n,] | lambda, phi);
    //target += -neg_binomial_2_lpmf(sum(counts[n,]) | alpha, phi_bar); //Note that alpha = sum(lambda)
    target += -lognormal_lpdf(sum(counts[n,]) | mean_log_bar, sd_log_bar);
  } 
}

// generated quantities {
//   vector[N] lnorm;
//   vector[N] nb;
//   
//   for(n in 1:N) {
//     lnorm[n] = -lognormal_lpdf(sum(counts[n,]) | mean_log_bar, sd_log_bar);
//     nb[n] = neg_binomial_2_lpmf(counts[n,]| lambda, phi);
//   } 
// }
