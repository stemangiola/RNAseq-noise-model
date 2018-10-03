data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure;

}

parameters {

  vector[G] lambda_z[N];
  vector[G] alpha;
  real<lower=0> sigma;

  real<lower=1> nu;
  real<lower=0> tau;

}
transformed parameters{
  vector[G] lambda[N];
  simplex[G] lambda_softmax[N];
  for(n in 1:N) lambda[n] = alpha + lambda_z[n] * sigma / sqrt(tau);
  for(n in 1:N) lambda_softmax[n] = softmax(lambda[n]);

}
model {

  // logNormal
  sigma ~ cauchy(0,2);
  tau ~ gamma(nu/2, nu/2);
  nu ~ gamma(2,0.1);

  sum(alpha) ~ normal(0,0.01 * G);
  alpha ~ normal(my_prior[1], my_prior[2]);

  // Multinomial
  for(n in 1:N) lambda_z[n] ~ normal(0,1);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(lambda_softmax[n]);

}
generated quantities{
  int<lower=0> counts_gen[N,G];

  for(n in 1:N) {
    counts_gen[n,] = multinomial_rng(lambda_softmax[n], exposure);
  }


}
