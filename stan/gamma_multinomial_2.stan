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
  int<lower=0, upper=1> generate_quantities;
  int<lower=0, upper=1> generate_log_lik;

  //Set to 1 for each sample that is held out
  int<lower=0, upper=1> holdout[N];

  int<lower=0> N_samples_log_lik;
}

transformed data {
  int<lower=0> N_gen = generate_quantities ? N : 0;
  int<lower=0> G_gen = generate_quantities ? G : 0;
  int<lower=0> N_log_lik = generate_log_lik ? N : 0;

  vector[2*G] Q_lambda = Q_sum_to_zero_QR(G);
  int<lower=1> exposure[N];

  for(n in 1:N) {
    exposure[n] = sum(counts[n,]);
  }
}


parameters {
  real<lower=0> lambda_sigma;
  vector[G - 1] lambda_raw;
  //real<lower=0> phi_raw_sigma;
  vector<lower=0>[G] theta[N];
  real<lower=0> asymptDisp;
  real<lower=0> extraPois;
}

transformed parameters {
  vector[G] lambda = sum_to_zero_QR(lambda_raw, Q_lambda) * lambda_sigma;
  vector<lower=0>[G] phi = asymptDisp + extraPois ./ exp(lambda);
  vector[G] gamma_rate = asymptDisp ./ exp(lambda) + extraPois;
}

model {
  for(n in 1:N) {
    if(holdout[n] == 0) {
      counts[n,] ~ multinomial(theta[n,] / sum(theta[n,]));
    }
    theta[n,] ~ gamma(phi, gamma_rate);
  }
  lambda_raw ~ normal(0, 1);
  lambda_sigma ~ normal(0, 2);
}

generated quantities{
  int<lower=0> counts_gen_naive[N_gen,G_gen];
  int<lower=0> counts_gen_geneWise[N_gen,G_gen];

  vector[G_gen] lambda_gen;
  vector[N_log_lik] log_lik;
  vector[N_log_lik] log_lik_sampled;

  if(generate_quantities) {
    // Sample gene wise rates
    //for(g in 1:G) lambda_gen[g] = normal_or_gammaLog_rng(lambda_mu, lambda_sigma, is_prior_asymetric);

    vector[G-1] lambda_gen_raw;
    matrix[N, G] theta_gen;

    for(g in 1:(G-1)) {
      lambda_gen_raw[g] = normal_rng(0, 1);
    }
    lambda_gen = sum_to_zero_QR(lambda_raw, Q_lambda) * lambda_sigma;

    {
      vector[G] gamma_rate_gen = asymptDisp ./ exp(lambda_gen) + extraPois;
      for(g in 1:G) {
        for(n in 1:N) {
          theta_gen[n, g] = gamma_rng(phi[g], gamma_rate_gen[g]);
        }
      }
    }

    // Sample gene wise sample wise abundances
    for(n in 1:N) {
      counts_gen_naive[n,] = multinomial_rng(to_vector(theta_gen[n,]) / sum(theta_gen[n,]), exposure[n]);
      counts_gen_geneWise[n,] = multinomial_rng(to_vector(theta[n,]) / sum(theta[n,]), exposure[n]);
    }

  }

    //log_lik for LOO
  if(generate_log_lik) {
    for(n in 1:N) {
      log_lik[n] = multinomial_lpmf(counts[n,] | to_vector(theta[n,]) / sum(theta[n,]));
    }

    {
      for(n in 1:N) {
        vector[N_samples_log_lik] log_lik_samp;
        for(s in 1:N_samples_log_lik) {
          vector[G] gamma_samp;
          for(g in 1:G) {
            gamma_samp[g] = gamma_rng(phi[g], gamma_rate[g]);
          }
          log_lik_samp[s] = multinomial_lpmf(counts[n,] | gamma_samp / sum(gamma_samp));
        }
        log_lik_sampled[n] = log_sum_exp(log_lik_samp) - log(N_samples_log_lik);
      }
    }
  }
}
