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

    real poisson_lindley_log_lpdf(int y, real mu, real alpha) {

      real exp_mu = exp(mu);
      real theta = 0.5/exp_mu[i]*(alpha+exp_mu[i]*(-1+sqrt((4*exp_mu[i]+square(alpha+exp_mu[i]))/square(exp_mu[i]))));
      real alpha_y = alpha + y;

      return (1 + alpha)*log(theta) - (1 + alpha_y)*log(1 + theta) + log(alpha + (alpha_y)/(1 + theta)) - lgamma(y+1) - lgamma(1 + alpha) + lgamma(alpha_y);
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
