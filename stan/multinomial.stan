data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real<lower=0> lambda_prior;
}

parameters {
  simplex[G] lambda_raw;
}

model {
  lambda_raw ~ dirichlet(rep_vector(lambda_prior, G));
  
  for(n in 1:N) {
    counts[n,] ~ multinomial(lambda_raw); 
  } 
}
