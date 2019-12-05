functions {
  vector lchoose_vec(vector n, real k) {
    vector[rows(n)] res;
    for(i in 1:rows(n)) {
      res[i] = lchoose(n[i],k);
    }
    return res;
  }

  real gwaring_lpmf(int[] y, real a, real k, real rho) {
    vector[size(y)] y_vec = to_vector(y);

    // return sum(lgamma(a + rho) + lgamma(k + rho)
    //           - lgamma(rho) //- lgamma(a + k + rho)
    //           + lgamma(a + y_vec) - lgamma(a)
    //           + lgamma(k + y_vec) - lgamma(k)
    //           - lgamma(a + k + rho + y_vec) //+ lgamma(a + k + rho)
    //           - lgamma(y_vec + 1)
    //           );
    return sum(
      -lchoose_vec(a + k + rho + y_vec - 1, a + rho - 1) - log(k + y_vec)
      + lchoose_vec(a + y_vec - 1, a - 1)
      + lchoose(k + rho - 1, k - 1) + log(rho)  );
  }

  real gwaring2_lpmf(int[] y, real mu, real phi, real r) {
      real s_border = sqrt((1 + phi) / (mu + phi));
      real s = s_border + inv(r);


      real a = phi + (1 + phi)/s;
      real k = phi + (mu + phi) * s;
      real rho = 2 + phi * (1 + s) * (1 + k) / (mu * s);

      return gwaring_lpmf(y | a, k, rho);
    }

  real gwaring4_lpmf(int[] y, real mu, real phi, real r) {
    real rho_low =  2 + phi + (phi * (1 + 2 * phi + 2 * sqrt((1 + phi) * (mu + phi)))) / mu;
    real rho = rho_low + (1 / r);

    if(rho > 1e06) {
       return neg_binomial_lpmf(y | mu, phi);
    } else {
      real kFirst = -phi + mu*(-2-phi+rho);
      real kSecond = sqrt(phi^2 + mu^2 * (2 + phi - rho)^2 - 2 * mu * phi * (-2 -3 * phi + rho + 2 * phi * rho));
      real k = (kFirst + kSecond) / (2 * phi);
      real a = mu * (rho - 1) / k;

      return gwaring_lpmf(y | a, k, rho);
    }
  }
}

data {
  int<lower=0> N;
  int y[N];
  real r_sd;
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
  y ~ gwaring4(mu, phi, r);

  mu ~ lognormal(3, 1);
  phi_raw ~ normal(0, 1);
  r ~ normal(0, r_sd);
  //r ~ lognormal(1, r_sd);
}

