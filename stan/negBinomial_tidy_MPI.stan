functions{

  vector lp_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {
   int M = xi[1];
   int N = xi[2];
   int S = xi[3];
   int D = xi[4];
   int G_per_shard = xi[5];
   int symbol_end[M+1] = xi[(5+1):(5+1+M)];
   int sample_idx[N] = xi[(5+1+M+1):(5+1+M+1+N-1)];
   int counts[N] = xi[(5+1+M+1+N):size(xi)];

   matrix[S,D] X = to_matrix(xr, S, D);

   matrix[D,G_per_shard] lambda_MPI = to_matrix(local_parameters[1:(G_per_shard*D)], G_per_shard, D)';
   vector[G_per_shard] sigma_MPI = local_parameters[((G_per_shard*D)+1):((G_per_shard*D) + G_per_shard)];
   vector[S] exposure_rate = local_parameters[((G_per_shard*D) + G_per_shard +1):((G_per_shard*D) + G_per_shard + S)];

   matrix[symbol_end[G_per_shard+1], G_per_shard] lambda_MPI_X = X[sample_idx[1:symbol_end[G_per_shard+1]]] * lambda_MPI;
   vector[G_per_shard] lp;

//print(lambda_MPI');

//print(sample_idx[1:symbol_end[G_per_shard+1]]);
  for(g in 1:G_per_shard){

   lp[g] =  neg_binomial_2_log_lpmf(

     counts[(symbol_end[g]+1):symbol_end[g+1]] |
     exposure_rate[sample_idx[(symbol_end[g]+1):symbol_end[g+1]]] +
     lambda_MPI_X[(symbol_end[g]+1):symbol_end[g+1],g],
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
	int<lower=0> D;
  int n_shards;
	int<lower=0> counts[n_shards, N];
	int<lower=0> symbol_end[n_shards, M+1];
	int<lower=0> sample_idx[n_shards, N];
	int<lower=0> G_per_shard[n_shards];
	int<lower=0> G_per_shard_idx[n_shards + 1];

	matrix<lower=0>[S,D] X;

}
transformed data {
  vector[0] global_parameters;
  real xr[n_shards, S * D];

  int<lower=0> int_MPI[n_shards, 5+(M+1)+N+N];

  // Shards - MPI
  for ( i in 1:n_shards ) xr[i,] = to_array_1d(X);
  for ( i in 1:n_shards ) {
  int M_N_Gps[5];
  M_N_Gps[1] = M;
  M_N_Gps[2] = N;
  M_N_Gps[3] = S;
  M_N_Gps[4] = D;
  M_N_Gps[5] = G_per_shard[i];

  int_MPI[i,] =
    append_array(
      append_array(
          append_array(M_N_Gps, symbol_end[i]),
        sample_idx[i]
      ),
      counts[i]
    );
  }
}

parameters {
  // Overall properties of the data
  real lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;
  real lambda_skew;

  // Gene-wise properties of the data
  matrix[G,D] lambda;
  vector[G] sigma_raw;

  vector[S] exposure_rate;

  // Signa linear model

  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;


}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ exp(sigma_raw);

// Shards - MPI
	vector[M * (D + 1)] lambda_sigma_MPI[n_shards];
	for( i in 1:(n_shards) ) {

	 // vector[ (M - G_per_shard[i]) * (D + 1) ] buffer = rep_vector(0.0,(M * (D + 1)) - (G_per_shard[i]*(D+1)));

		lambda_sigma_MPI[i] =
  		append_row(
  		  append_row(
  		    append_row(
        		to_vector(lambda[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])]),
      		  sigma[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])]
      		),
      		exposure_rate
      	),
      	rep_vector(0.0, (M - G_per_shard[i]) * (D + 1))
      );

	}

//	print(lambda);

}
model {

  // Overall properties of the data

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
  lambda[,1] ~  skew_normal(lambda_mu,lambda_sigma, lambda_skew);
  if(D>1) to_vector(lambda[,2:D]) ~ normal(0,1);
  sigma_raw ~ normal(sigma_slope * lambda[,1] + sigma_intercept,sigma_sigma);

	// Gene-wise properties of the data
	target += sum( map_rect( lp_reduce , global_parameters , lambda_sigma_MPI , xr , int_MPI ) );

}
generated quantities{
  vector[G] sigma_raw_gen;
  vector[G] lambda_gen;

  for(g in 1:G) sigma_raw_gen[g] = normal_rng(sigma_slope * lambda[g,1] + sigma_intercept,sigma_sigma);
  for(g in 1:G) lambda_gen[g] =  skew_normal_rng(lambda_mu,lambda_sigma, lambda_skew);


}
