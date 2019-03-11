functions{
	real gamma_log_lpdf(vector x_log, real a, real b){

		// This function is the  probability of the log gamma funnction
		// in case you have data that is aleady in log form

		vector[rows(x_log)] jacob = x_log; //jacobian
		real norm_constant = a * log(b) -lgamma(a);
		real a_minus_1 = a-1;
		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

	}


	vector gen_inv_logit(vector b1s, vector X_no_intercept, real inflection, real y_cross) {

		real b0 = -inflection * b1s[1];
		vector[num_elements(b1s) + 1] b = append_row(b0, b1s);
		matrix[rows(X_no_intercept), num_elements(b1s) + 1] X = append_col(rep_vector(1, rows(X_no_intercept)), X_no_intercept);

		return  y_cross * (1 + exp(-b0)) ./ to_vector(1 + exp(- (  X * b ))) ;
	}

}
data {
	int<lower=0> G;
	int<lower=0> S;
	int<lower=0> counts[G, S];


}
parameters {

	// Overall properties of the data
	real<lower=0> lambda_mu; // So is compatible with logGamma prior
	real<lower=0> lambda_sigma_raw;
	vector[S] exposure_rate;

	// Gene-wise properties of the data
	vector[G] lambda;
	vector<lower=0>[G] sigma_raw;

	// Signa linear model

	real sigma_inflection;
	real<lower=0> sigma_y_cross;
	vector<upper=0>[1] sigma_slope;
	real<lower=0>sigma_sigma;

}
transformed parameters {
	real<lower=0> lambda_sigma = lambda_sigma_raw / 1000;
	vector[G] sigma = 1.0 ./ sigma_raw;

	// Sigma linear model
	//vector[G] sigma_mu = exp(lambda * sigma_slope + sigma_0); //the expected values (linear predictor)
	vector[G] sigma_mu = gen_inv_logit(sigma_slope, lambda, sigma_inflection, sigma_y_cross) ;
	vector[G] sigma_alpha = sigma_mu .* sigma_mu / sigma_sigma; //shape parameter for the gamma distribution
	vector[G] sigma_beta = sigma_mu / sigma_sigma; //rate parameter for the gamma distribution

}
model {

	// Overall properties of the data
	lambda_mu ~ gamma(3,2);
	lambda_sigma_raw ~ normal(0,1);
	//sigma_raw ~ normal(0,1);
	exposure_rate ~ normal(0,1);
	sum(exposure_rate) ~ normal(0, 0.001 * S);

	sigma_inflection ~ normal(0,2);
	sigma_y_cross ~ cauchy(0,2);
	sigma_slope ~ normal(0,1);
	sigma_sigma ~ cauchy(0,2.5);

	// Gene-wise properties of the data
	for(g in 1:G) counts[g] ~ neg_binomial_2_log(lambda[g], sigma[g]);

  lambda ~ gamma_log_lpdf(lambda_mu, lambda_sigma);
	sigma_raw ~ gamma(sigma_alpha,sigma_beta);

	// Sample from data
	// if(omit_data==0) for(g in 1:G)
	// 	counts[symbol_start_end[g,1]:symbol_start_end[g,2]] ~ neg_binomial_2_log(
	// 		exposure_rate[sample_idx[symbol_start_end[g,1]:symbol_start_end[g,2]]] + lambda[g],
	// 		sigma[g]
	// 	);

}

