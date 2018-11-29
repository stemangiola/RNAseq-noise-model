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
	  else lpdf = normal_lpdf(x_log | 0, b);

	  return lpdf;
	}

	 real normal_or_gammaLog_rng(real a, real b, int is_prior_asymetric){
	  // This function takes care of the two prior choice
	  // without complicating too much the model itself
    real rng;
	  if(is_prior_asymetric == 1) rng = log(gamma_rng(a, b));
	  else rng = normal_rng(0, b);

	  return rng;
	}
}

data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure[N];

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;

  int<lower=0, upper=1> generate_quantities;

}

transformed data {
  int<lower=0> N_gen = generate_quantities ? N : 0;
  int<lower=0> G_gen = generate_quantities ? G : 0;
}

parameters {
  // Overall properties of the data
  real<lower=0> lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;

  // Gene-wise properties of the data
  vector[G] lambda;
}
model {

  // Overall properties of the data
  lambda_mu ~ normal(0,1);
  lambda_sigma ~ normal(0,2);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal_or_gammaLog(lambda_mu, lambda_sigma, is_prior_asymetric);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(lambda));

}
generated quantities{
  int<lower=0> counts_gen_naive[N_gen,G_gen];
  int<lower=0> counts_gen_geneWise[N_gen,G_gen];
  vector[N_gen] log_lik;

  vector[G_gen] lambda_gen;

  if(generate_quantities) {
    // Sample gene wise rates
    for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);

    // Sample gene wise sample wise abundances
    for(n in 1:N) {
      counts_gen_naive[n,] = multinomial_rng(softmax(lambda_gen), exposure[n]);
      counts_gen_geneWise[n,] = multinomial_rng(softmax(lambda), exposure[n]);
    }

    //log_lik for LOO
    for(n in 1:N) {
      log_lik[n] = multinomial_lpmf(counts[n,] | softmax(lambda));
    }
  }

}
