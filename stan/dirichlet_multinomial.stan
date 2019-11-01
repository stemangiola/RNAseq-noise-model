functions {

  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);

      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

	int[] dirichlet_multinomial_rng(vector alpha, int exposure) {
	    return multinomial_rng(dirichlet_rng(alpha), exposure);
	}

	real gamma_log_lpdf(vector x_log, real a){

    // This function is the  probability of the log gamma funnction
    // in case you have data that is aleady in log form

    // Constrain Mean to be 0
    real b = exp(digamma(a));

		vector[rows(x_log)] jacob = x_log; //jacobian
		real norm_constant = a * log(b) -lgamma(a);
		real a_minus_1 = a-1;
		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

	}

	real normal_or_gammaLog_lpdf(vector x_log, real a, int is_prior_asymetric){
	  // This function takes care of the two prior choice
	  // without complicating too much the model itself
    real lpdf;
	  if(is_prior_asymetric == 1) lpdf = gamma_log_lpdf(x_log | a);
	  else lpdf = normal_lpdf(x_log | 0, a);

	  return lpdf;
	}

	 real normal_or_gammaLog_rng(real a, int is_prior_asymetric){
	  // This function takes care of the two prior choice
	  // without complicating too much the model itself
    real rng;

    // Constrain Mean to be 0
    real b = exp(digamma(a));

	  if(is_prior_asymetric == 1) rng = log(gamma_rng(a, b));
	  else rng = normal_rng(0, a);

	  return rng;
	}
}

data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0> exposure[N];

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;

  //Set to 1 for each sample that is held out
  int<lower=0, upper=1> holdout[N];
}

parameters {

  // Overall properties of the data
  real<lower=0> lambda_prior; // So is compatible with logGamma prior
  real<lower=0> sigma;

  // Gene-wise properties of the data
  vector[G] lambda;
}
model {

  // Overall properties of the data
  lambda_prior ~ normal(0,1);
  sigma ~ gamma(3,2);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal_or_gammaLog(lambda_prior, is_prior_asymetric);

  // Sample from data
  for(n in 1:N) {
    if(holdout[n] == 0) {
      counts[n,] ~ dirichlet_multinomial(sigma * softmax(lambda));
    }
  }

}
generated quantities{
  int<lower=0> counts_gen_naive[N,G];
  int<lower=0> counts_gen_geneWise[N,G];
  vector[G] lambda_gen;
  vector[N] log_lik;

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_prior, is_prior_asymetric);

  // Sample gene wise sample wise abundances
  for(n in 1:N) {
    counts_gen_naive[n,] = dirichlet_multinomial_rng(sigma * softmax(lambda_gen), exposure[n]);
    counts_gen_geneWise[n,] = multinomial_rng(softmax(lambda), exposure[n]);
  }

  //log_lik for LOO
  for(n in 1:N) {
    log_lik[n] = dirichlet_multinomial_lpmf(counts[n,] | sigma * softmax(lambda));
  }

}
