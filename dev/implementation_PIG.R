
mu = 20
tau = 1


p0 = function() exp(tau^-1 * (1-(1+2*tau*mu)^(1/2)))
p1 = function(p_0) mu*(1+2*tau*mu)^-(1/2) * p(0)

p = function(y){

	switch(
		y %>% as.character,
		"0" = p0(),
		"1" = p1(),
		(2*tau*mu/(1+2*tau*mu)*(1-3/(2*y))*p(y-1)) + (mu^2/(1+2*tau*mu)*(1/(y*(y-1)))*p(y-2))
	)


}

p_optimised = function(y){

	p_arr =
		exp(tau^-1 * (1-(1+2*tau*mu)^(1/2)))

	if(y>0)
		p_arr =
			p_arr %>%
			c(
				mu*(1+2*tau*mu)^-(1/2) * .
			)

	if(y>1) foreach(y_dot = 2:y) %do% {
		p_arr = p_arr %>% c(
			(2*tau*mu/(1+2*tau*mu)*(1-3/(2*y_dot))* p_arr[y_dot]) + (mu^2/(1+2*tau*mu)*(1/(y_dot*(y_dot-1)))* p_arr[y_dot-1])
		)
	}

	p_arr %>% tail(n=1)

}

PIG = function(y, tau=1, mu=20){

	tau_mu_2 = tau*mu*2
	tau_mu_2_1 = tau_mu_2 + 1
	tau_mu_2_1_sqrt = tau_mu_2_1^(1/2)


	p_arr = exp(tau^-1 * (1-tau_mu_2_1_sqrt))

	if(y>0)	p_arr = p_arr %>% c( mu/tau_mu_2_1_sqrt * .	)

	if(y>1) foreach(y_dot = 2:y) %do% {
		y_dot_minus_1 = y_dot-1

		p_arr = p_arr %>% c(
			(tau_mu_2/(tau_mu_2_1)*(1-3/(2*y_dot))* p_arr[y_dot]) + (mu^2/tau_mu_2_1*(1/(y_dot*y_dot_minus_1))* p_arr[y_dot_minus_1])
		)
	}

	p_arr %>% tail(n=1) %>% log

}

PIG_log = function(y, tau=1, mu=20){

	tau_mu_2_log = log(tau) + log(mu) + log(2)
	tau_mu_2_1_log = log(exp(tau_mu_2_log) + 1)
	tau_mu_2_1_sqrt_log = 1/2 * tau_mu_2_1_log


	p_arr_log = tau^-1 * (1-exp(tau_mu_2_1_sqrt_log))

	if(y>0)	p_arr_log = p_arr_log %>% c( log(mu) - tau_mu_2_1_sqrt_log + . )

	if(y>1) foreach(y_dot = 2:y) %do% {
		y_dot_minus_1 = y_dot-1

		p_arr_log = p_arr_log %>% c(

			log(
				exp(  tau_mu_2_log - tau_mu_2_1_log + log(1-3/(2*y_dot)) + p_arr_log[y_dot] ) + exp( 2* log(mu) - tau_mu_2_1_log - log(y_dot) - log(y_dot_minus_1) + p_arr_log[y_dot_minus_1]  )
			)
		)
	}

	p_arr_log %>% tail(n=1)

}

PIG_log2 = function(y,  mu=20, tau=1){

	tau_mu_2_log = log(tau) + log(mu) + log(2)
	tau_mu_2_1_log = log(exp(tau_mu_2_log) + 1)
	tau_mu_2_1_sqrt_log = 1/2 * tau_mu_2_1_log


	p_arr_log = tau^-1 * (1-exp(tau_mu_2_1_sqrt_log))

	if(y>0)	p_arr_log = p_arr_log %>% c( log(mu) - tau_mu_2_1_sqrt_log + . )

	tau_mu_2_log_tau_mu_2_1_log = tau_mu_2_log - tau_mu_2_1_log
	two_log_mu_tau_mu_2_1_log = 2* log(mu) - tau_mu_2_1_log

	if(y>1) foreach(y_dot = 2:y) %do% {

		p_arr_log = p_arr_log %>% c(

			log(
				exp( tau_mu_2_log_tau_mu_2_1_log + log(1-3/(2*y_dot)) + p_arr_log[y_dot] ) +
				exp( two_log_mu_tau_mu_2_1_log - log(y_dot) - log(y_dot-1) + p_arr_log[y_dot-1]  )
			)
		)
	}

	p_arr_log %>% tail(n=1)

}

 poisson_inverse_gaussian_lpmf = function( y,  log_mu,  tau){


   tau_mu_2_log = log(tau) + log_mu + 0.6931472; # log(2)
   tau_mu_2_1_log =  log(exp(tau_mu_2_log) + 1);
   tau_mu_2_1_sqrt_log = 1.0/2 * tau_mu_2_1_log;
 # vector[y + 1] p_arr_log;
  # Start to create array
   p_arr_log = array()
  p_arr_log[1] = 1/tau * (1-exp(tau_mu_2_1_sqrt_log));

  if(y>0)	p_arr_log[2] = log_mu - tau_mu_2_1_sqrt_log + p_arr_log[1];

  if(y>1) {
     tau_mu_2_log_tau_mu_2_1_log = tau_mu_2_log - tau_mu_2_1_log;
     two_log_mu_tau_mu_2_1_log = 2* log_mu - tau_mu_2_1_log;

    for(y_dot in 2:y) {



      p_arr_log[y_dot + 1] =
        log(
          exp(tau_mu_2_log_tau_mu_2_1_log + log(1-3.0/(2*y_dot)) + p_arr_log[y_dot]) +
          exp(two_log_mu_tau_mu_2_1_log - log(y_dot) - log(y_dot-1) + p_arr_log[y_dot-1])
        );
    }
  }

   p_arr_log[ y + 1 ];

}

mu = 3
tau = 2

PIG_log2(10, exp(3), 2)
poisson_inverse_gaussian_lpmf(10, 3, 2)

system.time(PIG_log2(1000))

p_short = function(y){
	log(1/factorial(y)) + log(p0()) + ifelse(y>0, foreach(y_dot = 0:(y-1), .combine=sum) %do% log((y_dot +1) * ) , 0)
}


p(0)
p(1)
p(2)
p(3)
p(4)
p(5)
p(20)

system.time( p(20) )


sapply(0:1000, function(s) system.time( p_optimised(s) )[3])
system.time( p_optimised2(20) )[3]
