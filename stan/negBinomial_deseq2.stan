data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> counts[N,G];

  //1 - logmean, 2- logsd of the counts after normalization
  real my_prior[2];

  //Set to 1 for each sample that is held out
  int<lower=0, upper=1> holdout[N];

  int<lower=0, upper=1> generate_quantities;
  int<lower=0, upper=1> generate_log_lik;

  real normalization[N];
}

transformed data {
  int<lower=0> N_gen = generate_quantities ? N : 0;
  int<lower=0> G_gen = generate_quantities ? G : 0;
  int<lower=0> N_log_lik = generate_log_lik ? N : 0;
}


parameters {
  // Overall properties of the data

  // Gene-wise log-means
  vector[G] lambda;
  real<lower=0> asymptDisp;
  real<lower=0> extraPois;
}
transformed parameters {
  vector<lower=0>[G] phi = asymptDisp + extraPois ./ exp(lambda);
}

model {
  lambda ~ normal(my_prior[1], my_prior[2]);


  // Sample from data
  for(n in 1:N) {
    if(holdout[n] == 0) {
      counts[n,] ~ neg_binomial_2_log(lambda + log(normalization[n]), phi);
    }
  }

}
generated quantities{
  int<lower=0> counts_gen_geneWise[N_gen,G_gen];
  vector[N_log_lik] log_lik;


  if(generate_quantities) {
    // Sample gene wise sample wise abundances
    for(n in 1:N) {
      for(g in 1:G) {
        counts_gen_geneWise[n,g] = neg_binomial_2_log_rng(lambda[g] + log(normalization[n]), phi[g]);
      }
    }
  }

  if(generate_log_lik) {
    for(n in 1:N) {
      log_lik[n] = neg_binomial_2_log_lpmf(counts[n,] | lambda + log(normalization[n]), phi);
    }
  }

}
