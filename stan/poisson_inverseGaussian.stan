functions{

real poisson_inverse_gaussian_log_lpmf(int[] y, real log_mu, real sigma){

  int y_max = max(y);
	real sigma_mu_2_log = log(sigma) + log_mu + 0.6931472; // log(2)
	real sigma_mu_2_1_log =  log1p_exp(sigma_mu_2_log);
	real sigma_mu_2_1_sqrt_log = 1.0/2 * sigma_mu_2_1_log;
	vector[y_max + 1] p_arr_log;
	int y_plus_1[size(y)];

  for(i in 1:size(y)) y_plus_1[i] = y[i] +1;

	// Start to create array
	p_arr_log[1] = 1.0/sigma * (1-exp(sigma_mu_2_1_sqrt_log));

	if(y_max>0)	p_arr_log[2] = log_mu - sigma_mu_2_1_sqrt_log + p_arr_log[1];
	if(y_max>1) {
	  real sigma_mu_2_log_sigma_mu_2_1_log = sigma_mu_2_log - sigma_mu_2_1_log;
	  real two_log_mu_sigma_mu_2_1_log = 2* log_mu - sigma_mu_2_1_log;

	  for(y_dot in 2:y_max) {

  		p_arr_log[y_dot + 1] =
  			log_sum_exp(
  				 sigma_mu_2_log_sigma_mu_2_1_log + log(1-3.0/(2*y_dot)) + p_arr_log[y_dot],
  				 two_log_mu_sigma_mu_2_1_log - log(y_dot) - log(y_dot-1) + p_arr_log[y_dot-1]
  			);
	  }
	}

  return sum(p_arr_log[ y_plus_1 ]);

}

    // int poisson_inverse_gamma_log_rng(int y, real mu, real beta, real sigma){
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
    //   return poisson(inv_gamma_rng(.., ..)); // real mu, real beta, real sigma
    //
    // }
}
data {
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0> sigma;
  real mu;
}
model {
mu ~ normal(0,5);
sigma ~ normal(0,1);
y ~ poisson_inverse_gaussian_log(mu,sigma);

// print( poisson_inverse_gaussian_lpmf(x | log(100), 1));

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
