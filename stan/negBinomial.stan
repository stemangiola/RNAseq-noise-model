data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
}

parameters {

  // Overall properties of the data
  real lambda_mu;
  real<lower=0> lambda_sigma;
  real<lower=0> sigma_raw;
  real exposure_rate[N];

  // Gene-wise properties of the data
  vector[G] lambda;

}
transformed parameters {
  real sigma = 1/sqrt(sigma_raw);
}
model {

  // Overall properties of the data
  lambda_mu ~ normal(0,5);
  lambda_sigma ~ normal(0,2);
  sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.01 * N);

  // Gene-wise properties of the data
  lambda ~ normal(lambda_mu, lambda_sigma);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ neg_binomial_2_log(exposure_rate[n] + lambda, sigma);

}
generated quantities{
  int<lower=0> counts_gen_naive[N,G];
  int<lower=0> counts_gen_geneWise[N,G];
  vector[G] lambda_gen;

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_rng(lambda_mu, lambda_sigma);

  // Sample gene wise sample wise abundances
  for(n in 1:N) for(g in 1:G) {
    counts_gen_naive[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda_gen[g], sigma);
    counts_gen_geneWise[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda[g],  sigma);
  }


}
