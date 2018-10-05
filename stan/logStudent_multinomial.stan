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
  real<lower=0> sigma;
  real<lower=1> nu;
  real<lower=0> tau;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector[G] theta_z[N];

}
transformed parameters{
  vector[G] theta[N];
  for(n in 1:N) theta[n] = lambda + theta_z[n] * sigma / sqrt(tau);
}
model {

  // Overall properties of the data
  lambda_mu ~ normal(0,1);
  lambda_sigma ~ cauchy(0,2);
  sigma ~ cauchy(0,2);
  tau ~ gamma(nu/2, nu/2);
  nu ~ gamma(2,0.1);

  // Gene-wise properties of the data
  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal(lambda_mu, lambda_sigma);
  for(n in 1:N) theta_z[n] ~ normal(0,1);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(theta[n]));

}
generated quantities{
  int<lower=0> counts_gen[N,G];
  vector[G] lambda_gen;
  vector[G] theta_gen[N];

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_rng(lambda_mu, lambda_sigma);
  for(n in 1:N) for(g in 1:G) theta_gen[n,g] = student_t_rng(nu, lambda_gen[g],sigma);

  // Sample gene wise sample wise abundances
  for(n in 1:N) {
    counts_gen[n,] = multinomial_rng(softmax(theta_gen[n]), exposure[n]);
  }

}


