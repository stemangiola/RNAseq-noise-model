functions{
  	real gamma_log_lpdf(real x_log, real a, real b){

      // This function is the  probability of the log gamma funnction
      // in case you have data that is aleady in log form

  		real jacob = x_log; //jacobian
  		real norm_constant = a * log(b) -lgamma(a);
  		real a_minus_1 = a-1;
  		return  jacob  + norm_constant * 1 + (  x_log * a_minus_1 - exp(x_log) * b ) ;

  	}

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
  ordered[2] lambda_mu; // So is compatible with logGamma prior
  vector<lower=0>[2] lambda_sigma;
  real<lower=0, upper=1> lambda_prop;

  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector<lower=0>[G] sigma_raw;

  // Signa linear model

  real sigma_inflection;
  real<lower=0> sigma_y_cross;
  vector<upper=0>[1] sigma_slope;
  real<lower=0>sigma_sigma;

}
transformed parameters {

  //real<lower=0> lambda_beta = lambda_beta_raw / 1000;
  vector[G] sigma = 1.0 ./ sigma_raw;

  // Sigma linear model
  //vector[G] sigma_mu = exp(lambda * sigma_slope + sigma_0); //the expected values (linear predictor)
  vector[G] sigma_mu = gen_inv_logit(sigma_slope, lambda, sigma_inflection, sigma_y_cross) ;
  vector[G] sigma_alpha = sigma_mu .* sigma_mu / sigma_sigma; //shape parameter for the gamma distribution
  vector[G] sigma_beta = sigma_mu / sigma_sigma; //rate parameter for the gamma distribution

  vector[2] lambda_alpha = lambda_mu .* lambda_mu ./ lambda_sigma / 10;
  vector[2] lambda_beta = lambda_mu ./ lambda_sigma / 10;


}
model {

  // Overall properties of the data
  lambda_mu[1] ~ cauchy(0,2.5);
  lambda_mu[2] ~ normal(4,2);
  lambda_sigma ~ normal(0,1);
  lambda_prop ~ beta(1, 10);

  //sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  sigma_inflection ~ normal(0,2);
  sigma_y_cross ~ cauchy(0,2);
  sigma_slope ~ normal(0,1);
  sigma_sigma ~ cauchy(0,2.5);

  // Gene-wise properties of the data
  for(g in 1:G)
    target +=
      log_mix(
        lambda_prop,
        gamma_log_lpdf(lambda[g] | lambda_alpha[1], lambda_beta[1]),
        gamma_log_lpdf(lambda[g] | lambda_alpha[2], lambda_beta[2])
      );

  sigma_raw ~ gamma(sigma_alpha,sigma_beta);

  // Sample from data
  if(omit_data==0) for(g in 1:G)
      counts[symbol_start_end[g,1]:symbol_start_end[g,2]] ~ neg_binomial_2_log(
        exposure_rate[sample_idx[symbol_start_end[g,1]:symbol_start_end[g,2]]] + lambda[g],
        sigma[g]
      );

}
generated quantities{
  vector<lower=0>[G] sigma_raw_gen;

  for(g in 1:G) sigma_raw_gen[g] = gamma_rng(sigma_alpha[g],sigma_beta[g]);


}
// generated quantities{
//   int<lower=0> counts_gen_naive[S,G];
//   int<lower=0> counts_gen_geneWise[S,G];
//   vector[G] lambda_gen;
//
//   // Sample gene wise rates
//   for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_alpha, lambda_beta, is_prior_asymetric);
//
//   // Sample gene wise sample wise abundances
//   for(s in 1:S) for(g in 1:G) {
//     counts_gen_naive[s,g] = neg_binomial_2_log_rng(exposure_rate[s] + lambda_gen[g], sigma[g]);
//     counts_gen_geneWise[s,g] = neg_binomial_2_log_rng(exposure_rate[s] + lambda[g],  sigma[g]);
//   }
//
//
// }
