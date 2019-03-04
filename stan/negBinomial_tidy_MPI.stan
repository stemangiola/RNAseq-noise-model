functions{

	 vector lp_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {
	 	int M = xi[1];
	 	int N = xi[2];
	 	int S = xi[3];
	 	int G_per_shard = xi[4];
	 	int symbol_start[M+1] = xi[(4+1):(4+1+M)];
	 	int sample_idx[N] = xi[(4+1+M+1):(4+1+M+1+N-1)];
	 	int counts[N] = xi[(4+1+M+1+N):size(xi)];


	 	vector[G_per_shard] lambda_MPI = local_parameters[1:G_per_shard];
	 	vector[G_per_shard] sigma_MPI = local_parameters[(G_per_shard+1):(G_per_shard*2)];
	 	vector[S] exposure_rate = local_parameters[((M*2)+1):rows(local_parameters)];

	 	vector[G_per_shard] lp;


	 for(g in 1:G_per_shard){
	 	lp[g] =  neg_binomial_2_log_lpmf(
	 	  counts[symbol_start[g]:symbol_start[g+1]-1] |
	 	  exposure_rate[sample_idx[symbol_start[g]:symbol_start[g+1]-1]] +
	 	  lambda_MPI[g],
	 	  sigma_MPI[g]
	 	 );
	 }


    return [sum(lp)]';
  }

}
data {
  int<lower=0> N;
  int<lower=0> M;
	int<lower=0> G;
	int<lower=0> S;
  int n_shards;
	int<lower=0> counts[n_shards, N];
	int<lower=0> symbol_start[n_shards, M+1];
	int<lower=0> sample_idx[n_shards, N];
	int<lower=0> G_per_shard[n_shards];
	int<lower=0> G_per_shard_idx[n_shards + 1];

}
transformed data {
  vector[0] global_parameters;
  real xr[n_shards, 0];

  int<lower=0> int_MPI[n_shards, 4+(M+1)+N+N];

  // Shards - MPI
  for ( i in 1:n_shards ) {
  int M_N_Gps[4];
  M_N_Gps[1] = M;
  M_N_Gps[2] = N;
  M_N_Gps[3] = S;
  M_N_Gps[4] = G_per_shard[i];

  int_MPI[i,] = append_array(append_array(append_array(M_N_Gps, symbol_start[i]), sample_idx[i]), counts[i]);

  }
}
parameters {
  // Overall properties of the data
  real lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;
  real lambda_skew;
  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G] lambda;
  vector<lower=0>[G] sigma_raw;

  // Signa linear model

  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;

}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ sigma_raw;

// Shards - MPI
	vector[2*M] lambda_sigma_MPI[n_shards];
	for( i in 1:(n_shards) ) {

	  vector[ (M*2) - (G_per_shard[i]*2) ] buffer = rep_vector(0.0,(M*2) - (G_per_shard[i]*2));

		lambda_sigma_MPI[i] =
  		append_row(
  		  append_row(
  		    append_row(
      		  lambda[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])],
      		  sigma[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])]
      		),
      		buffer
      	),
      	exposure_rate
      );

	}

}
model {

  // Overall properties of the data
int i = 1;
  lambda_mu ~ normal(0,2);
  lambda_sigma ~ normal(0,2);
  lambda_skew ~ normal(0,2);

  //sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

  // Gene-wise properties of the data
  // lambda ~ normal_or_gammaLog(lambda_mu, lambda_sigma, is_prior_asymetric);
  lambda ~  skew_normal(lambda_mu,lambda_sigma, lambda_skew);
  sigma_raw ~ lognormal(sigma_slope * lambda + sigma_intercept,sigma_sigma);

	// Gene-wise properties of the data
	target += sum( map_rect( lp_reduce , global_parameters , lambda_sigma_MPI , xr , int_MPI ) );

}
generated quantities{
  vector<lower=0>[G] sigma_raw_gen;
  vector[G] lambda_gen;

  for(g in 1:G) sigma_raw_gen[g] = lognormal_rng(sigma_slope * lambda[g] + sigma_intercept,sigma_sigma);
  for(g in 1:G) lambda_gen[g] =  skew_normal_rng(lambda_mu,lambda_sigma, lambda_skew);


}
