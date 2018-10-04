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
}
model {

  lambda ~ normal(my_prior[1], my_prior[2]);
  sigma ~ gamma(2,3);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ neg_binomial_2_log(lambda, sigma);

}
generated quantities{
  int<lower=0> counts_gen[N,G];

  for(n in 1:N) for(g in 1:G) {
    counts_gen[n,g] = neg_binomial_2_log_rng(lambda[g], sigma);
  }


}
