functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);

      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
	 	int[] dirichlet_multinomial_rng(vector alpha, int exposure) {
	    return multinomial_rng(dirichlet_rng(alpha), exposure);
	  }

}




data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0> exposure;

}

parameters {

  vector[G] lambda;
  real<lower=1> sigma;
}
transformed parameters{
  simplex[G] lambda_softmax = softmax(lambda);
}
model {

  sum(lambda) ~ normal(0,0.01 * G);
  lambda ~ normal(my_prior[1], my_prior[2]);
  lambda ~ gamma(3, 2);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(sigma * lambda_softmax);

}
generated quantities{
  int<lower=0> counts_gen[N,G];

  for(n in 1:N) {
    counts_gen[n,] = dirichlet_multinomial_rng(lambda_softmax, exposure);
  }

}
