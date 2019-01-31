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
  int<lower=0> N;
  int<lower=0> G;
  matrix[N,G] direct_proportions;
  real<lower=0> lambda_prior;
  //real<lower=0, upper=1> softmax_vals[N,G];
}

transformed data {
  vector[2*N] Q_r = Q_sum_to_zero_QR(N);
}


parameters {

  // Overall properties of the data
  real<lower=0> sigma;

  // Gene-wise properties of the data
  vector[G] lambda;
  matrix[N - 1,G] theta_z_raw;

}
transformed parameters{
  matrix[N, G] theta_z;
  matrix[N, G] theta;
  for(g in 1:G) theta_z[,g] = sum_to_zero_QR(theta_z_raw[,g], Q_r);
  for(n in 1:N) theta[n] = to_row_vector(lambda) + theta_z[n] * sigma;
}

model {

  // Overall properties of the data
  lambda_prior ~ normal(0,1);
  sigma ~ normal(0,1);

  // Gene-wise properties of the data
  //sum(lambda) ~ normal(0,0.001 * G);
  lambda ~ normal(0, lambda_prior);

  //for(g in 1:G) sum(theta_z[,g]) ~ normal(0, 0.001 * N);
  //for(n in 1:N) sum(theta[n]) ~ normal(0, 0.001 * G);
  //for(n in 1:N) sum(theta_z[n]) ~ normal(0, 0.001 * G);

  for(n in 1:(N - 1)) {
    theta_z_raw[n] ~ normal(0,1);
  }

  for(n in 1:N) direct_proportions[n,] ~ normal(theta[n,],0.1);
  // Sample from data
  //if(omit_data==0) for(n in 1:N) counts[n,] ~ multinomial(softmax(theta[n]));

}

