functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	real alpha_plus = sum(alpha);

      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
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
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];
  real<lower=0> lambda_raw_prior;

  //0 - Implemented as DM, 1 - implemented as NB, 2 - implemented as gamma_multinom test
  int<lower=0,upper=2> variant;

  int<lower=0,upper=1> theta_given;

  vector<lower=0>[theta_given ? 1 : 0] theta_data;
  vector<lower=0>[theta_given ? 0 : 1] theta_prior_log_mean;
  vector<lower=0>[theta_given ? 0 : 1] theta_prior_log_sd;

  int<lower=0,upper=1> alpha_given;
  vector<lower=0>[alpha_given ? 1 : 0] alpha_data;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_mean;
  vector<lower=0>[alpha_given ? 0 : 1] alpha_prior_log_sd;

}

transformed data {
  int<lower=0> sums[N];
  real dummy_x_r[0];

  for(n in 1:N) {
    sums[n] = sum(counts[n,]);
  }
}

parameters {
  simplex[G] lambda_raw;
  vector<lower=0>[theta_given ? 0 : 1] log_theta_param_raw;
  vector<lower=0>[alpha_given ? 0 : 1] log_alpha_param_raw;
}

transformed parameters {
  vector<lower=0>[theta_given ? 0 : 1] theta_param = exp(log_theta_param_raw .* theta_prior_log_sd + theta_prior_log_mean);
  vector<lower=0>[alpha_given ? 0 : 1] alpha_param = exp(log_alpha_param_raw .* alpha_prior_log_sd + alpha_prior_log_mean);
  real alpha = alpha_given ? alpha_data[1] : alpha_param[1];
  vector[G] lambda = alpha * lambda_raw;
  real theta = theta_given ? theta_data[1] : theta_param[1];
}

model {
  //lambda_raw ~ dirichlet(rep_vector(lambda_raw_prior, G));
  lambda_raw ~ lognormal(0, lambda_raw_prior);

  if(!theta_given) {
    log_theta_param_raw ~ normal(0,1);
  }
  if(!alpha_given) {
    log_alpha_param_raw ~ normal(0,1);
  }

  for(n in 1:N) {
    if(variant == 0) {
      counts[n,] ~ dirichlet_multinomial_lpmf(lambda * sums[n] / theta);
    } else if(variant == 1) {
      counts[n,] ~ neg_binomial_2(lambda * sums[n], lambda * sums[n] / theta);
      target += -neg_binomial_lpmf(sum(counts[n,]) | alpha * sums[n], alpha * sums[n] / theta);//note that alpha = sum(lambda)
    } else if(variant == 2) {
      counts[n,] ~ gamma_multinomial_nb(log(lambda * sums[n]), lambda * sums[n] /theta, dummy_x_r);
    } else {
      reject("Invalid variant");
    }
  }
}
