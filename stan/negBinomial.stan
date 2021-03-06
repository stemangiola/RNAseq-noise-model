functions{
  	real gamma_log_lpdf(vector x_log, real a, real b){

      // This function is the  probability of the log gamma funnction
      // in case you have data that is aleady in log form

  		vector[rows(x_log)] jacob = x_log; //jacobian
  		real norm_constant = a * log(b) -lgamma(a);
  		real a_minus_1 = a-1;
  		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

  	}

  	real normal_or_gammaLog_lpdf(vector x_log, real a, real b, int is_prior_asymetric){
  	  // This function takes care of the two prior choice
  	  // without complicating too much the model itself
      real lpdf;
  	  if(is_prior_asymetric == 1) lpdf = gamma_log_lpdf(x_log | a, b);
  	  else lpdf = normal_lpdf(x_log | a, b);

  	  return lpdf;
  	}

  	real normal_or_gammaLog_rng(real a, real b, int is_prior_asymetric){
  	  // This function takes care of the two prior choice
  	  // without complicating too much the model itself
      real rng;
  	  if(is_prior_asymetric == 1) rng = log(gamma_rng(a, b));
  	  else rng = normal_rng(a, b);

  	  return rng;
  	}
}
data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  //Set to 1 for each sample that is held out
  int<lower=0, upper=1> holdout[N];

  int<lower=0, upper=1> generate_quantities;
  int<lower=0, upper=1> generate_log_lik;

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;
}

transformed data {
  int<lower=0> N_gen = generate_quantities ? N : 0;
  int<lower=0> G_gen = generate_quantities ? G : 0;
  int<lower=0> N_log_lik = generate_log_lik ? N : 0;
}


parameters {

  // Overall properties of the data
  real<lower=0> lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma_raw;
  real<lower=0> sigma_raw;
  real exposure_rate[N];

  // Gene-wise properties of the data
  vector[G] lambda;

}
transformed parameters {
  real<lower=0> lambda_sigma = lambda_sigma_raw / 1000;
  real sigma = 1/sqrt(sigma_raw);
}
model {

  // Overall properties of the data
  lambda_mu ~ gamma(3,2);
  lambda_sigma_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * N);

  // Gene-wise properties of the data
  lambda ~ normal_or_gammaLog(lambda_mu, lambda_sigma, is_prior_asymetric);

  // Sample from data
  for(n in 1:N) {
    if(holdout[n] == 0) {
      counts[n,] ~ neg_binomial_2_log(exposure_rate[n] + lambda, sigma);
    }
  }

}
generated quantities{
  int<lower=0> counts_gen_naive[N_gen,G_gen];
  int<lower=0> counts_gen_geneWise[N_gen,G_gen];
  vector[N_log_lik] log_lik;

  vector[G_gen] lambda_gen;

  if(generate_quantities) {
    // Sample gene wise rates
    for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);

    // Sample gene wise sample wise abundances
    for(n in 1:N) for(g in 1:G) {
      counts_gen_naive[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda_gen[g], sigma);
      counts_gen_geneWise[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda[g],  sigma);
    }
  }

  if(generate_log_lik) {
    for(n in 1:N) {
      log_lik[n] = neg_binomial_2_log_lpmf(counts[n,] | exposure_rate[n] + lambda, sigma);
    }
  }

}
