data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure[N];

}

parameters {

 // Overall properties of the data
  real<lower=0> lambda_sigma;
  real<lower=0> sigma;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector[G] theta_z[N];

}
transformed parameters{
  vector[G] theta[N];
  for(n in 1:N) theta[n] = lambda + theta_z[n] * sigma;
}
model {

  // Overall properties of the data
  lambda_sigma ~ normal(0,2);
  sigma ~ cauchy(0,2);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal(0, lambda_sigma);
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
  for(g in 1:G) lambda_gen[g] = normal_rng(0, lambda_sigma);
  for(n in 1:N) for(g in 1:G) theta_gen_naive[n,g] = normal_rng(lambda_gen[g],sigma);
  for(n in 1:N) for(g in 1:G) theta_gen_geneWise[n,g] = normal_rng(lambda[g],sigma);

  // Sample gene wise sample wise abundances
  for(n in 1:N) {
    counts_gen_naive[n,] = multinomial_rng(softmax(theta_gen_naive[n]), exposure[n]);
    counts_gen_geneWise[n,] = multinomial_rng(softmax(theta_gen_geneWise[n]), exposure[n]);
  }

}
