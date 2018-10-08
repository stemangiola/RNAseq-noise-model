
functions {
// Generate a vector of size N that sums to zero and marginal
// distributions for all elements are N(0,1)
// Original code due to andre.pfeuffer
// https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884/31

   vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;
    real scale = 1/sqrt(1-inv(N));

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0)) * scale;
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1)) * scale;
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
}

data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real my_prior[2];
  int<lower=0, upper=1> omit_data;
  int<lower=0> exposure[N];
  int<lower=0, upper=1> sigma_prior_type;
  real<lower=0> sigma_prior_sigma;
}

transformed data {
  vector[2*G] Q_r = Q_sum_to_zero_QR(G);
}

parameters {
  // Overall properties of the data
  real lambda_mu;
  real<lower=0> lambda_sigma;

  // Gene-wise properties of the data
  vector[G - 1] lambda_raw;
}

transformed parameters {
  vector[G] lambda = sum_to_zero_QR(lambda_raw, Q_r) * lambda_sigma;
}

model {

  // Overall properties of the data
  lambda_mu ~ normal(0,5);
  if(sigma_prior_type == 0) {
    lambda_sigma ~ cauchy(0, sigma_prior_sigma);
  } else {
    lambda_sigma ~ normal(0, sigma_prior_sigma);
  }

  lambda_raw ~ normal(0,1);

  // Sample from data
  if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(lambda));

}
generated quantities{
  int<lower=0> counts_gen_naive[N,G];
  int<lower=0> counts_gen_geneWise[N,G];

  vector[G] lambda_gen;

  // Sample gene wise rates
  for(g in 1:G) lambda_gen[g] = normal_rng(lambda_mu, lambda_sigma);

  // Sample gene wise sample wise abundances
  for(n in 1:N) {
    counts_gen_naive[n,] = multinomial_rng(softmax(lambda_gen), exposure[n]);
    counts_gen_geneWise[n,] = multinomial_rng(softmax(lambda), exposure[n]);
  }


}
