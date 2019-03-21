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

  real max_s(vector mus, vector phis) {
    return min(append_row(-log(mus) + log(mus + phis), log(phis ./ mus + 1)));
  }

  real s_transform(real s, vector mus, vector phis) {
    return -exp(s) + max_s(mus, phis);
  }

  vector nb_sum_log_Kd_eq(vector y, vector theta, real[] x_r, int[] x_i) {
    int G = rows(theta) / 2;
    vector[G] log_mus = theta[1:G];
    vector[G] mus = exp(log_mus);
    vector[G] phis = theta[(G + 1) : (2 * G)];
    real s = s_transform(y[1], mus, phis);
    real sum_y = x_i[1];
    real value = log_sum_exp(log(phis) + log_mus + s - log(phis + mus - mus * exp(s))) - log(sum_y);
    //print("sum_y = ",sum_y,"; log_mus = c(", log_mus, "); phis = c(", phis, "); s = ", s, "; value = ", value);
    return to_vector({value});
  }

  real neg_binomial_sum_lpmf(int sum_y, vector log_mus, vector phis, real[] dummy_x_r) {
    int G = rows(log_mus);
    vector[G] mus = exp(log_mus);

    //print(sum_y,", ", log_mus, ", ", phis);
    //{

    // Solve the saddlepoint equation
    vector[2 * G] solver_params = append_row(log_mus, phis);
    vector[1] solver_guess = to_vector({0});
    vector[1] s_vec_raw = algebra_solver(nb_sum_log_Kd_eq, solver_guess, solver_params, dummy_x_r,  {sum_y});
    real s = s_transform(s_vec_raw[1], mus, phis);
    //Calculate the saddlepoint mass
    vector[G] log_denominator_s = log(phis + mus - mus * exp(s));
    real K_s = sum(phis .* (log(phis) - log_denominator_s));
    real log_Kdd_s = log_sum_exp(log(phis) + log_mus + log(phis + mus) + s - 2 * log_denominator_s);
    real sum_lpmf = -0.5 * (log(2*pi()) + log_Kdd_s) + K_s - s * sum_y;

    return sum_lpmf;
    //}
  }

  real gamma_multinomial_nb_lpmf(int[] y, vector log_mus, vector phis, real[] dummy_x_r) {
    // real sum_mus = sum(mus);
    // real phi_bar = (sum_mus .* sum_mus) ./ sum((mus .* mus) ./ phis);
    // vector[size(y)] y_phis = to_vector(y) + phis;
    // vector[size(y)] yp1 = to_vector(y) + 1;
    //int sum_y = sum(y);
    //real log_sum_y = log(sum_y);

    return neg_binomial_2_log_lpmf(y | log_mus, phis) - neg_binomial_sum_lpmf(sum(y)| log_mus, phis, dummy_x_r);
    // return lgamma(sum_y + 1) - lgamma(sum_y + phi_bar) + sum(lgamma(y_phis) - lgamma(yp1)) +
    //   (sum_y + phi_bar) * log(sum_mus + phi_bar) - sum_y * log(sum_mus) +
    //   sum(to_vector(y) .* log(mus) - (y_phis) .* log(mus + phis)) +
    //   lgamma(phi_bar) - sum(lgamma(phis)) +
    //   sum(phis .* log(phis)) - phi_bar * log(phi_bar);
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
}

transformed data {
  int<lower=0> N_gen = generate_quantities ? N : 0;
  int<lower=0> G_gen = generate_quantities ? G : 0;
  int<lower=0> N_log_lik = generate_log_lik ? N : 0;

  vector[2*G] Q_lambda = Q_sum_to_zero_QR(G);
  int<lower=1> exposure[N];
  real dummy_x_r[0];

  for(n in 1:N) {
    exposure[n] = sum(counts[n,]);
  }
}


parameters {
  real<lower=0> lambda_sigma;
  vector[G - 1] lambda_raw;
  //real<lower=0> phi_raw_sigma;
  real<lower=0> asymptDisp;
  real<lower=0> extraPois;
}

transformed parameters {
  vector[G] lambda = sum_to_zero_QR(lambda_raw, Q_lambda) * lambda_sigma;
  vector[G] exp_lambda = exp(lambda);
  vector<lower=0>[G] phi = asymptDisp + extraPois ./ exp_lambda;
}

model {
  for(n in 1:N) {
    if(holdout[n] == 0) {
      counts[n,] ~ gamma_multinomial_nb_lpmf(lambda, phi, dummy_x_r);
    }
  }
  lambda_raw ~ normal(0, 1);
  lambda_sigma ~ normal(0, 2);
}

generated quantities{
  int<lower=0> counts_gen_naive[N_gen,G_gen];
  int<lower=0> counts_gen_geneWise[N_gen,G_gen];

  vector[G_gen] lambda_gen;
  vector[N_log_lik] log_lik;

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
      //counts_gen_geneWise[n,] = multinomial_rng(to_vector(theta[n,]) / sum(theta[n,]), exposure[n]);
    }

  }

    //log_lik for LOO
  if(generate_log_lik) {
    for(n in 1:N) {
      log_lik[n] = gamma_multinomial_nb_lpmf(counts[n,] | lambda, phi, dummy_x_r);
    }
  }


}
