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
}
transformed parameters{
  simplex[G] lambda_softmax = softmax(lambda);
}
model {

  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal(my_prior[1], my_prior[2]);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(lambda_softmax);

}
generated quantities{
  int<lower=0> counts_gen[N,G];

  for(n in 1:N) {
    counts_gen[n,] = multinomial_rng(lambda_softmax, exposure);
  }


}
