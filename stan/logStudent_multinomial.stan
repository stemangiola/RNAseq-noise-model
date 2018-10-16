functions{
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
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure[N];

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;
  int<lower=1> variance_size;

}

parameters {

  // Overall properties of the data
  real<lower=0> lambda_prior; // So is compatible with logGamma prior
  real<lower=1> nu;
  real<lower=0> tau;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector<lower=0>[variance_size] sigma_raw;
  vector[G] theta_z[N];

}
transformed parameters{
  vector<lower=0>[G] sigma = variance_size == 1 ? rep_vector(sigma_raw[1], G) : sigma_raw;
  vector[G] theta[N];
  for(n in 1:N) theta[n] = lambda + theta_z[n] .* sigma / sqrt(tau);
}
model {

  // Overall properties of the data
  lambda_prior ~ normal(0,1);
	sigma_raw ~ normal(0,1);
  tau ~ gamma(nu/2, nu/2);
  nu ~ gamma(2,0.1);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.001 * G);
  lambda ~ normal_or_gammaLog(lambda_prior, is_prior_asymetric);
  for(n in 1:N) sum(theta_z[n]) ~ normal(0, 0.001 * G);
  for(n in 1:N) theta_z[n] ~ normal(0,1);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(theta[n]));

}
generated quantities{
  int<lower=0> counts_gen_naive[N,G];
  int<lower=0> counts_gen_geneWise[N,G];
  vector[G] lambda_gen;
  vector[G] theta_gen_naive[N];
  vector[G] theta_gen_geneWise[N];

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_prior, is_prior_asymetric);
  for(n in 1:N) for(g in 1:G) theta_gen_naive[n,g] = student_t_rng(tau, lambda_gen[g],sigma[g]);
  for(n in 1:N) for(g in 1:G) theta_gen_geneWise[n,g] = student_t_rng(tau, lambda[g],sigma[g]);

  // Sample gene wise sample wise abundances
  for(n in 1:N) counts_gen_naive[n,] = multinomial_rng(softmax(theta_gen_naive[n]), exposure[n]);
  for(n in 1:N) counts_gen_geneWise[n,] = multinomial_rng(softmax(theta_gen_geneWise[n]), exposure[n]);
}
