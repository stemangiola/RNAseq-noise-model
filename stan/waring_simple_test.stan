functions {
  real waring2_lpmf(int[] y, real mu, real phi, real r) {
      real s_border = 0.5 * (sqrt(4*mu + 4 * mu * phi + phi^2) / mu - phi / mu);
      real tau = (s_border + r) * mu / phi;

      real k = (tau + 1)*phi;
      real rho = (phi + 1) / tau + phi + 2;
      real a = mu .* (rho - 1) ./ k;
      vector[size(y)] y_vec = to_vector(y);

      return sum(lgamma(a + rho) + lgamma(k + rho)
                - lgamma(rho) - lgamma(a + k + rho)
                + lgamma(a + y_vec) - lgamma(a)
                + lgamma(k + y_vec) - lgamma(k)
                - lgamma(a + k + rho + y_vec) + lgamma(a + k + rho)
                - lgamma(y_vec + 1)
                );

    }
}

data {
  int<lower=0> N;
  int y[N];
}

parameters {
  real<lower=0> mu;
  real<lower=0> phi_raw;
  real<lower=0> r;
}

transformed parameters {
  real phi = 1/sqrt(phi_raw);
}

model {
  y ~ waring2(mu, phi, r);

  mu ~ lognormal(3, 1);
  phi_raw ~ normal(0, 1);
  r ~ normal(0, 5);
}

