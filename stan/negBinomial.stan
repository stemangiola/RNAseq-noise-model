data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure;

}

parameters {

  vector[G] lambda;
  real<lower=1> sigma;
  real<lower=0> sigma_raw;
  real exposure_rate[N];
}
transformed parameters {
  real sigma = 1/sqrt(sigma_raw);
}
model {

  lambda ~ normal(my_prior[1], my_prior[2]);
  sigma ~ gamma(2,3);
  sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.01 * N);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ neg_binomial_2_log(exposure_rate[n] + lambda, sigma);

}
generated quantities{
  int<lower=0> counts_gen[N,G];

  for(n in 1:N) for(g in 1:G) {
    counts_gen[n,g] = neg_binomial_2_log_rng(lambda[g], sigma);
    counts_gen[n,g] = neg_binomial_2_log_rng(exposure_rate[n] + lambda_gen[g], sigma);
  }


}
