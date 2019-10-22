functions{

  real poisson_tweedie_raw_lpmf(int[] x, real a, real b, real c);

  real poisson_tweedie_lpmf (int[] x, real log_mu, real D, real a) {
    real b = (exp(log_mu) * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a));
    real c = (D - 1)/(D - a);
    return poisson_tweedie_raw_lpmf(x | a, b, c);
  }

  real poisson_tweedie_abc_lpmf (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw_lpmf(x | a, b, c);
  }
}

data{
  int<lower=0> N;
  int y[N];
  real<lower=0> mu;
  real<lower=0> D;
}
parameters{
  real<lower=0> a_raw;
}

transformed parameters {
  real<upper=1> a = 1 - a_raw;
}

model {
  a_raw ~ lognormal(1, 1);
  y ~ poisson_tweedie(log(mu), D, a);
}

