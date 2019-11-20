functions{
  	real gamma_log_lpdf(vector x_log, real a, real b){

      // This function is the  probability of the log gamma function
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

  	vector gen_inv_logit(vector b1s, vector X_no_intercept, real inflection, real y_cross) {

  	  real b0 = -inflection * b1s[1];
  	  vector[num_elements(b1s) + 1] b = append_row(b0, b1s);
  	  matrix[rows(X_no_intercept), num_elements(b1s) + 1] X = append_col(rep_vector(1, rows(X_no_intercept)), X_no_intercept);

      return  y_cross * (1 + exp(-b0)) ./ to_vector(1 + exp(- (  X * b ))) ;
    }


}
data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> S;
  int<lower=0> counts[N];
  int<lower=0> sample_idx[N];

  int<lower=0> symbol_start_end[G,2];
  int<lower=0, upper=1> omit_data;

  // Alternative models
  int<lower=0, upper=1> is_prior_asymetric;
}
parameters {

  // Overall properties of the data
  real lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;
  real lambda_skew;
  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector<lower=0>[G] sigma_raw;

  // Signa linear model

  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;

}
transformed parameters {
  vector[G] sigma = 1.0 ./ sigma_raw;

  // Sigma linear model
  //vector[G] sigma_mu = exp(lambda * sigma_slope + sigma_0); //the expected values (linear predictor)
  // vector[G] sigma_mu = gen_inv_logit(sigma_slope, lambda, sigma_inflection, sigma_y_cross) ;
  // vector[G] sigma_alpha = sigma_mu .* sigma_mu / sigma_sigma / 10; //shape parameter for the gamma distribution
  // vector[G] sigma_beta = sigma_mu / sigma_sigma / 10; //rate parameter for the gamma distribution

}
model {

  // Overall properties of the data
  lambda_mu ~ normal(0,2);
  lambda_sigma ~ normal(0,2);
  lambda_skew ~ normal(0,2);

  //sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

  // Gene-wise properties of the data
  // lambda ~ normal_or_gammaLog(lambda_mu, lambda_sigma, is_prior_asymetric);
  lambda ~  skew_normal(lambda_mu,lambda_sigma, lambda_skew);
  sigma_raw ~ lognormal(sigma_slope * lambda + sigma_intercept,sigma_sigma);

  // Sample from data
  if(omit_data==0) for(g in 1:G)
      counts[symbol_start_end[g,1]:symbol_start_end[g,2]] ~ neg_binomial_2_log(
        exposure_rate[sample_idx[symbol_start_end[g,1]:symbol_start_end[g,2]]] + lambda[g],
        sigma[g]
      );

}
generated quantities{
  vector<lower=0>[G] sigma_raw_gen;

  for(g in 1:G) sigma_raw_gen[g] = lognormal_rng(sigma_slope * lambda[g] + sigma_intercept,sigma_sigma);


}
// generated quantities{
//   int<lower=0> counts_gen_naive[S,G];
//   int<lower=0> counts_gen_geneWise[S,G];
//   vector[G] lambda_gen;
//
//   // Sample gene wise rates
//   for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);
//
//   // Sample gene wise sample wise abundances
//   for(s in 1:S) for(g in 1:G) {
//     counts_gen_naive[s,g] = neg_binomial_2_log_rng(exposure_rate[s] + lambda_gen[g], sigma[g]);
//     counts_gen_geneWise[s,g] = neg_binomial_2_log_rng(exposure_rate[s] + lambda[g],  sigma[g]);
//   }
//
//
// }
