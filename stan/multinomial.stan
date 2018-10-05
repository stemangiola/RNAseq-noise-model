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
  real lambda_mu;
  real<lower=0> lambda_sigma;

  // Gene-wise properties of the data
  vector[G] lambda;
}
model {

  // Overall properties of the data
  lambda_mu ~ normal(0,1);
  lambda_sigma ~ cauchy(0,2);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal(lambda_mu, lambda_sigma);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(lambda));

}
generated quantities{
  int<lower=0> counts_gen[N,G];
  vector[G] lambda_gen;

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_rng(lambda_mu, lambda_sigma);

  // Sample gene wise sample wise abundances
  for(n in 1:N) {
    counts_gen[n,] = multinomial_rng(softmax(lambda_gen), exposure[n]);
  }


}
