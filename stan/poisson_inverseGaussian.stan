functions{
  // 	real gamma_log_lpdf(vector x_log, real a, real b){
  //
  //     // This function is the  probability of the log gamma funnction
  //     // in case you have data that is aleady in log form
  //
  // 		vector[rows(x_log)] jacob = x_log; //jacobian
  // 		real norm_constant = a * log(b) -lgamma(a);
  // 		real a_minus_1 = a-1;
  // 		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;
  //
  // 	}
  //
  // 	real normal_or_gammaLog_lpdf(vector x_log, real a, real b, int is_prior_asymetric){
  // 	  // This function takes care of the two prior choice
  // 	  // without complicating too much the model itself
  //     real lpdf;
  // 	  if(is_prior_asymetric == 1) lpdf = gamma_log_lpdf(x_log | a, b);
  // 	  else lpdf = normal_lpdf(x_log | a, b);
  //
  // 	  return lpdf;
  // 	}
  //
  // 	real normal_or_gammaLog_rng(real a, real b, int is_prior_asymetric){
  // 	  // This function takes care of the two prior choice
  // 	  // without complicating too much the model itself
  //     real rng;
  // 	  if(is_prior_asymetric == 1) rng = log(gamma_rng(a, b));
  // 	  else rng = normal_rng(a, b);
  //
  // 	  return rng;
  // 	}

real poisson_inverse_gaussian_lpmf(int y, real log_mu, real tau){


	real tau_mu_2_log = log(tau) + log_mu + 0.6931472; // log(2)
	real tau_mu_2_1_log =  log(exp(tau_mu_2_log) + 1);
	real tau_mu_2_1_sqrt_log = 1.0/2 * tau_mu_2_1_log;
	vector[y + 1] p_arr_log;

	// Start to create array
	p_arr_log[1] = 1.0/tau * (1-exp(tau_mu_2_1_sqrt_log));

	if(y>0)	p_arr_log[2] = log_mu - tau_mu_2_1_sqrt_log + p_arr_log[1];

	if(y>1) {
	  real tau_mu_2_log_tau_mu_2_1_log = tau_mu_2_log - tau_mu_2_1_log;
	  real two_log_mu_tau_mu_2_1_log = 2* log_mu - tau_mu_2_1_log;

	  for(y_dot in 2:y) {



  		p_arr_log[y_dot + 1] =
  			log_sum_exp(
  				 tau_mu_2_log_tau_mu_2_1_log + log(1-3.0/(2*y_dot)) + p_arr_log[y_dot],
  				 two_log_mu_tau_mu_2_1_log - log(y_dot) - log(y_dot-1) + p_arr_log[y_dot-1]
  			);
	  }
	}

	return p_arr_log[ y + 1 ];

}

    // int poisson_inverse_gamma_log_rng(int y, real mu, real beta, real tau){
    //
    //   // Sample inverse gaussian
    //   real y = normal_rng(0,1)^2;
    //   real mu_2 = mu^2
    //   real x = mu + (mu_2*y)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*y + mu_2*y*y)
    //   real z = uniform(0,1);
    //
    //   if (z <= (mu)/(mu + x)) return x;
    //   else return mu_2/x;
    //
    //   return poisson(inv_gamma_rng(.., ..)); // real mu, real beta, real tau
    //
    // }
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
transformed data{
    print(poisson_inverse_gaussian_lpmf(10 | 3, 2));
}
parameters {

  // Overall properties of the data
  real<lower=0> lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma_raw;
  //real<lower=0> sigma;
  real<lower=0, upper=1> tau;
  //real exposure_rate[N];

  // Gene-wise properties of the data
  vector[G] lambda_raw;


}
transformed parameters {
  real<lower=0> lambda_sigma = lambda_sigma_raw;
  //if(is_prior_asymetric) lambda_sigma = lambda_sigma / 1000;
  vector[G] lambda = lambda_mu + lambda_raw * lambda_sigma;



}
model {

  // Overall properties of the data
  lambda_mu ~ gamma(3,2);
  lambda_sigma_raw ~  gamma(6,5);
  //sigma ~ normal(0,2);
  tau ~ normal(0,2);
  //exposure_rate ~ normal(0,1);
  //sum(exposure_rate) ~ normal(0, 0.001 * N);

  // Gene-wise properties of the data
  lambda_raw ~ normal(0, 1);

  // Sample from data
   if(omit_data==0) for(n in 1:N) for(g in 1:G) counts[n,g] ~  poisson_inverse_gaussian(lambda[g], tau);

}
// generated quantities{
//   int<lower=0> counts_gen_naive[N,G];
//   int<lower=0> counts_gen_geneWise[N,G];
//   vector[G] lambda_gen;
//
//   // Sample gene wise rates
//   for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);
//
//   // Sample gene wise sample wise abundances
//   for(n in 1:N) for(g in 1:G) {
//     counts_gen_naive[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda_gen[g], sigma);
//     counts_gen_geneWise[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda[g],  sigma);
//   }
//
//
// }
