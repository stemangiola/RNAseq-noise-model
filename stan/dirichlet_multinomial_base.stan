functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);

      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}


data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real<lower=0> lambda_raw_prior;

  //0 - Implemented as DM, 1 - implemented as NB
  int<lower=0,upper=1> variant;

  int<lower=0,upper=1> theta_given;

  vector<lower=0>[theta_given ? 1 : 0] theta_data;
  vector<lower=0>[theta_given ? 0 : 1] theta_prior_log_mean;
  vector<lower=0>[theta_given ? 0 : 1] theta_prior_log_sd;

  int<lower=0,upper=1> alpha_given;
  vector<lower=0>[alpha_given ? 1 : 0] alpha_data;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_mean;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_sd;

}

transformed data {
  int<lower=0> sums[N];

  for(n in 1:N) {
    sums[n] = sum(counts[n,]);
  }
}

parameters {
  simplex[G] lambda_raw;
  vector<lower=0>[theta_given ? 0 : 1] log_theta_param_raw;
  vector<lower=0>[alpha_given ? 0 : 1] log_alpha_param_raw;
}

transformed parameters {
  vector<lower=0>[theta_given ? 0 : 1] theta_param = exp(log_theta_param_raw .* theta_prior_log_sd + theta_prior_log_mean);
  vector<lower=0>[alpha_given ? 0 : 1] alpha_param = exp(log_alpha_param_raw .* alpha_prior_log_sd + alpha_prior_log_mean);
  real alpha = alpha_given ? alpha_data[1] : alpha_param[1];
  vector[G] lambda = alpha * lambda_raw;
  real theta = theta_given ? theta_data[1] : theta_param[1];
}

model {
  //lambda_raw ~ dirichlet(rep_vector(lambda_raw_prior, G));
  lambda_raw ~ lognormal(0, lambda_raw_prior);

  if(!theta_given) {
    log_theta_param_raw ~ normal(0,1);
  }
  if(!alpha_given) {
    log_alpha_param_raw ~ normal(0,1);
  }

  for(n in 1:N) {
    if(variant == 0) {
      counts[n,] ~ dirichlet_multinomial_lpmf(lambda * sums[n] / theta);
    } else {
      counts[n,] ~ neg_binomial_2(lambda, lambda / theta);
      target += -neg_binomial_lpmf(sum(counts[n,]) | alpha * sums[n], alpha * sums[n] / theta);//note that alpha = sum(lambda)
    }
  }
}
