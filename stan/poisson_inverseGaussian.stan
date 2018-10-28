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

  	real poisson_inverse_gaussian_log_original_lpdf(int y, real mu, real beta, real tau) {

  	  real exp_mu = exp(mu);
  	  real exp_mu_beta_tau = exp_mu*beta+tau;
  	  real exp_mu_beta_tau_exp_mu_2 = exp_mu_beta_tau*exp_mu*2.0;
  	  real prob;

      if(y==0) prob =
        exp(1/exp_mu_beta_tau*(1-(1+exp_mu_beta_tau_exp_mu_2)^0.5));

      if(y==1) prob =
        exp_mu*(1+exp_mu_beta_tau_exp_mu_2)^-0.5 *
        exp(poisson_inverse_gamma_log_lpdf(0,mu,beta, tau));

      prob =
        exp_mu_beta_tau_exp_mu_2 /
        (1+exp_mu_beta_tau_exp_mu_2) *
        (1. - 3. /(2.*y)) *
        exp(poisson_inverse_gamma_log_lpdf(y-1,mu,beta, tau)) +
        exp_mu^2 /
        (1+exp_mu_beta_tau_exp_mu_2) /
        (y*(y-1)) *
        exp(poisson_inverse_gamma_log_lpdf(y-2,mu,beta, tau));

      return log(prob);
    }

  	real poisson_inverse_gamma_log_lpdf(int y, real mu, real beta, real tau) {

  	  real mu = mu;
  	  real mu_beta_tau = mu + log(beta+tau);
  	  real mu_beta_tau_mu_2 = mu_beta_tau + mu + 2.0;
  	  real prob;

      if(y==0) prob =
        1 / mu_beta_tau * (1-(1+mu_beta_tau_mu_2)^0.5);

      if(y==1) prob =
        mu -
        ( 0.5 * log1p(mu_beta_tau_mu_2)) +
        poisson_inverse_gamma_log_lpdf(0,mu,beta, tau);

      prob =
      log_sum_exp(
        log(mu_beta_tau_mu_2) -
        log1p(mu_beta_tau_mu_2) +
        log_sum_exp(0.0, - ( log(3.) - log(2.) + log(y) ) ) +
        poisson_inverse_gamma_log_lpdf(y-1,mu,beta, tau),

        ( 2 * mu ) -
        log1p(mu_beta_tau_mu_2) -
        log(y) + log(y-1) +
        poisson_inverse_gamma_log_lpdf(y-2,mu,beta, tau)
      );

      return prob;
    }

    int poisson_inverse_gamma_log_rng(int y, real mu, real beta, real tau){

      return poisson(inv_gamma_rng(.., ..)); // real mu, real beta, real tau

    }
}
data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;
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
  if(omit_data==0) for(n in 1:N) counts[n,] ~ poisson_inverse_gaussian_log(exposure_rate[n] + lambda, sigma, tau);

}
generated quantities{
  int<lower=0> counts_gen_naive[N,G];
  int<lower=0> counts_gen_geneWise[N,G];
  vector[G] lambda_gen;

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);

  // Sample gene wise sample wise abundances
  for(n in 1:N) for(g in 1:G) {
    counts_gen_naive[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda_gen[g], sigma);
    counts_gen_geneWise[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda[g],  sigma);
  }


}
