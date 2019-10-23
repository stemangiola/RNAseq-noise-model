functions{

  real poisson_tweedie_raw(int[] x, real a, real b, real c);
  real poisson_tweedie_raw_da(int[] x, real a, real b, real c);


  real poisson_tweedie_lpmf (int[] x, real mu, real D, real a) {
    real b = (mu * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a));
    real c = (D - 1)/(D - a);
    return poisson_tweedie_raw(x, a, b, c);
  }

  real poisson_tweedie_abc_lpmf (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw(x, a, b, c);
  }

  real poisson_tweedie_abc_wrapper (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw(x, a, b, c);
  }

  real poisson_tweedie_abc_da_wrapper (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw_da(x, a, b, c);
  }

}

data{
  int<lower=0> N;
  int y[N];
  real mu_prior_logmean;
  real<lower=0> mu_prior_sigma;
  real Dm1_prior_logmean;
  real<lower=0> Dm1_prior_sigma;
  //real<lower=0> b;
  //real<lower=0> c;
}
parameters{
  real<lower=0, upper = 1> a;
  real<lower=0> mu;
  real<lower=1> D;
}

transformed parameters {
  //real<upper=1> a = 1 - a_raw;
}

model {
  //a_raw ~ lognormal(1, 1);
  mu ~ lognormal(mu_prior_logmean, mu_prior_sigma);
  (D - 1) ~ lognormal(Dm1_prior_logmean, Dm1_prior_sigma);
  y ~ poisson_tweedie(mu, D, a);
  //target += poisson_tweedie_raw(y, a, b, c);
}

