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
      real x_aux_tmp = x_aux;
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux_tmp + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
}

data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=0> counts[N, G];
  vector<lower=0>[G] alpha;
}

transformed data {
  vector[2*G] Q_lambda = Q_sum_to_zero_QR(G);
}


parameters {
  vector[G - 1] lambda_raw;
//  vector[G] lambda;
  real<lower=0> phi_raw;
  vector<lower=0>[G] theta[N];
}

transformed parameters {
  real<lower=0> phi = 1/sqrt(phi_raw);
  vector[G] lambda = sum_to_zero_QR(lambda_raw, Q_lambda);
}

model {
  for(n in 1:N) {
    counts[n,] ~ multinomial(theta[n,] / sum(theta[n,]));
  }
  for(g in 1:G) {
    theta[, g] ~ gamma(phi, phi / (alpha[g] * exp(lambda[g])));
  }
  phi_raw ~ normal(0, 1);
  lambda_raw ~ normal(0, 1);

  // lambda ~ normal(0, 1);
  // sum(lambda) ~ normal(0, 0.01 * G);
}
