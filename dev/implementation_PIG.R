
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

PIG_alt = function(y, alpha=1, mu=20){


  w = (mu^2 + alpha^2)^(1/2) - mu
  gamma = -1/2

  p_arr = exp(w - alpha)

  if(y>0)	p_arr = p_arr %>% c(mu * w / alpha * .	)

  if(y>1) foreach(y_dot = 2:y) %do% {
    y_dot_minus_1 = y_dot-1

    p_arr = p_arr %>% c(

      ( 2* mu * w / alpha^2 * (y + gamma - 1) / y * p_arr[y_dot] ) +
      ( (mu*w/alpha)^2 / (y*(y-1)) * p_arr[y_dot_minus_1] )

    )
  }

  p_arr %>% tail(n=1) %>% log

}

PIG_log = function(y,  mu=20, tau=1){

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



approximated_modified_bessel_second_kind_log = function(z, v, s = max(0,v-10)){

  0.5 * ( log(pi) - log(2) ) + (-z - 0.5 * log(z)) +


  foreach(j = s:floor(v-0.5), .combine = c) %do% {

    lgamma(j + abs(v) - 0.5 + 1) -

    ( lgamma(j + 1) + lgamma(-j + abs(v) - 0.5 + 1) ) -

    j*(log(2) + log(z))

  } %>% matrixStats::logSumExp()
}

approximated_further_modified_bessel_second_kind_log = function(z, v, s = max(0,v-10)){
  foreach(j = s:floor(v-0.5), .combine = c) %do% {
    ( (j + abs(v) - 0.5 + 1) %>% lgamma() ) -
      ( lgamma(j + 1) + lgamma(-j + abs(v) - 0.5 + 1) ) -
      j*(log(2) + log(z))
  } %>% c(


  ) %>% matrixStats::logSumExp()
}

approximated_modified_bessel_second_kind_log2 = function(z,v){
  v = v - 1/2
  1/2 * (log(pi) - log(2) - log(v)) + -v * log(exp(1)*z / (2*v))
}




mu = 3
tau = 2

PIG_log(10, exp(3), 2)
PIG_log2(10, exp(3), 2)
poisson_inverse_gaussian_lpmf(10, 3, 2)

system.time(PIG_log2(1000))




# Bessel function
library(doParallel)
registerDoParallel()
foreach(alpha = 1:30,  .combine = bind_rows) %dopar% {
  foreach(v = 1:500,  .combine = bind_rows) %do% {
    tibble(
      v = v,
      alpha = alpha,
      error =
        approximated_modified_bessel_second_kind_log(alpha, v, 0) -
        approximated_modified_bessel_second_kind_log(alpha, v)
    )
  }
} %>%
  {
    (.) %>%
    ggplot(aes(x = alpha, y = v)) +
    geom_tile(aes(fill = error), color = "white") +
    coord_equal() +
    # geom_text(
    #   aes(label = error, 2),
    #   size = 4,
    #   colour = "black"
    # ) +
    scale_fill_gradient(
      low = "#c3e87d",
      high = "#f29496"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

    (.)
  } %>%
  {
    (.) %>%
      filter(alpha == 30) %>%
      ggplot(aes(x=v, y=error)) +
      geom_point()
  }

foreach(alpha = 1000,  .combine = bind_rows) %do% {
  foreach(j = 1:1000,  .combine = bind_rows) %dopar% {
    tibble(
      j = j,
      alpha = alpha,
      alpha_term = j*(log(2) + log(alpha)),
      v_term =     ( (j + abs(v) - 0.5 + 1) %>% lgamma() ) -
        ( lgamma(j + 1) + lgamma(-j + abs(v) - 0.5 + 1) ) ,
      # tot_term = ( (j + abs(v) - 0.5 + 1) %>% lgamma() ) -
      #   ( lgamma(j + 1) + lgamma(-j + abs(v) - 0.5 + 1) ) -
      #   j*(log(2) + log(alpha)),
      full = approximated_modified_bessel_second_kind_log(alpha, j, 0),
      approx = approximated_modified_bessel_second_kind_log(alpha, j)
    )
  }
} %>%
gather(term, value, -j, -alpha) %>%
group_by(term) %>%
arrange(j) %>%
mutate(value_cum = log(cumsum(exp(value)))) %>%
ungroup() %>%
{
  (.) %>%
    ggplot(aes(x=j, y=value, color=term)) +
    geom_line()
}

approximated_modified_bessel_second_kind_log2(0.5, 1000)

system.time(approximated_modified_bessel_second_kind_log(0.5, 1000, 0) )
system.time(approximated_modified_bessel_second_kind_log(0.5, 1000, 990) )
system.time(approximated_modified_bessel_second_kind_log2(0.5, 1000) )


my_ppoisinvgauss = function(x, mu, phi){

  sqrt(1/phi) * sqrt(1/(pi/2)) * exp(1/(phi * mu))/factorial(x) *
    (sqrt(2 *phi *(1 + (2 * phi * mu^2)^(-1))))^(-(x - 0.5)) *
    besselK( sqrt(2/phi *(1 + (2*  phi * mu^2)^(-1))), x - 0.5)

  #
  # sqrt(2/(pi * dispersion))*exp(1/(dispersion * mean))/factorial(x) *
  #   sqrt(2 * dispersion * (1 + 1/(2 * dispersion * mean^2)))^(-(x - 1/2)) *
  #   besselK(x - 1/2, sqrt(2/dispersion * (1 + 1/(2*  dispersion  * mean^2))))
}

my_ppoisinvgauss_log = function(x, mu, phi){

   0.5 * -log(phi)  + 0.5 * (- log(pi) + log(2)) + 1/(phi * mu) - lgamma(x + 1 ) +


    (-(x - 0.5)) * ( 0.5 * (log(2) + log(phi) + log1p(1/ (2 * phi * mu^2))) )  +

   #log( (sqrt(2 *phi *(1 + (2 * phi * mu^2)^(-1))))^(-(x - 0.5)) ) +

   log( besselK( sqrt(2/phi *(1 + (2*  phi * mu^2)^(-1))), x - 0.5) )

  #
  # sqrt(2/(pi * dispersion))*exp(1/(dispersion * mean))/factorial(x) *
  #   sqrt(2 * dispersion * (1 + 1/(2 * dispersion * mean^2)))^(-(x - 1/2)) *
  #   besselK(x - 1/2, sqrt(2/dispersion * (1 + 1/(2*  dispersion  * mean^2))))
}


my_ppoisinvgauss(100, 100,  1/1000) %>% log
my_ppoisinvgauss_log(100, 100,  1/1000)
actuar::dpoisinvgauss(100, mean = 100, dispersion = 1/1000) %>% log

sapply(0:100, actuar::dpoisinvgauss, 100, 1000) %>% plot

approximated_modified_bessel_second_kind_log(0.1, 100, 0)
besselK(0.1, 100) %>% log
approximated_modified_bessel_second_kind_log(10, 100, 0)
besselK(10, 100) %>% log
approximated_modified_bessel_second_kind_log(30, 100, 0)
besselK(30, 100) %>% log
approximated_modified_bessel_second_kind_log(30, 1000, 0)
besselK(30, 1000) %>% log



 my_dpoisinvgauss = function( x,  mu,  phi)
{

     phim = phi * mu
     lphi = log(phi);
     a = 1/(2 * phim * mu)
     y = x - 0.5;
     M_LN_SQRT_PId2 = 0.225791352644727432363097614947
     M_LN2 = 0.693147180559945309417232121458

    logA = -lphi/2 - M_LN_SQRT_PId2 + 1/phim;
    logB = (M_LN2 + lphi + log1p(a))/2;
    lgamma1p = function(x) lgamma( x + 1)

  lpx = logA - y * logB - lgamma1p(x);
 # K = besselK(exp(logB - lphi), y);
  K_log = approximated_modified_bessel_second_kind_log(exp(logB - lphi), y, 0)
  # lpx + log(K);
  lpx + K_log
 }

 actuar::dpoisinvgauss(100, 100, 0.01, log =  T)
my_dpoisinvgauss(500, 100, 1/2000)

tbl =
  seq(0,1000, 10) %>%
  tibble(
    x = .,
    `dpoisinvgauss` = sapply(. , actuar::dpoisinvgauss, 100, 2000, log=T),
    `dpoisinvgauss + stable besselK` = sapply(. , my_dpoisinvgauss, 100, 1/2000),
    `dnbinom` = sapply(. , dnbinom,mu= 100, size=15)%>% log,
    `dSICHEL` = sapply(. , gamlss.dist::dSICHEL,mu= 100, sigma=0.1, nu=-30, log=T),
    dnorm = sapply(. ,  metRology::dt.scaled, df = 40, mean=100, sd=50) %>% log ,
    `dt.scaled` = sapply(x , dnorm,mean= 100, 50)%>% log
  ) %>%
  gather(model, `log density`, -x)

tbl  %>%
  dplyr::filter(model %in% c("dpoisinvgauss", "dpoisinvgauss + stable besselK")) %>%
  ggplot(aes(x = x, y=`log density`, color = model)) + geom_point() + facet_grid(~model) + my_theme +
  ggtitle("sapply(seq(0,1000, 10) , actuar::dpoisinvgauss, 100, 2000, log=T)")

tbl  %>%
  dplyr::filter(model != "dpoisinvgauss") %>%
  dplyr::mutate(
    Type = ifelse(
      model %in% c("dpoisinvgauss", "dpoisinvgauss + stable besselK", "dnbinom", "dSICHEL"),
      "Discrete",
      "Continuous"
    )
  ) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Discrete", "Continuous"))) %>%
  dplyr::filter(Type == "Discrete") %>%
  ggplot(aes(x = x, y=`log density`, color = model)) + geom_point() + facet_grid(~ Type) + my_theme +
  ggtitle("NB vs. PIG vs. SICHEL")



 my_tofySICHEL2 = function(y, mu, sigma, nu, lbes, cvec, ans, ny, maxy) {
  maxyp1 = maxy + 1;
  tofY = c();

  for (i in 1:ny) {
    iy = y[i]+1;
    tofY[1] = (mu[i]/cvec[i])*((1+2*sigma[i]*mu[i]/cvec[i])^-0.5)*exp(lbes[i]);
    alpha = sqrt(1 + 2*sigma[i]*mu[i]/cvec[i])/sigma[i];
    sumT = 0;
    foreach (j = 2:iy) %do% {
      browser()
      tofY[j] = ( cvec[i]*sigma[i]*(2*(j+nu[i])/mu[i]) + (1/tofY[j-1])) * ((mu[i]/(sigma[i]*alpha*cvec[i])^2) );
      print(tofY[j])
      sumT = sumT + log(tofY[j-1]);
    }
    ans[i] = sumT;
  }

  ans
  }





my_dSICHEL = function (x, mu = 1, sigma = 1, nu = -0.5, log = FALSE)
{

  ly <- max(length(x), length(mu), length(sigma), length(nu))
  x <- rep(x, length = ly)
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)
  nu <- rep(nu, length = ly)
  cvec <- exp(log(besselK((1/sigma), nu + 1)) - log(besselK((1/sigma), nu)))
  alpha <- sqrt(1 + 2 * sigma * (mu/cvec))/sigma
  lbes <- log(besselK(alpha, nu + 1)) - log(besselK((alpha),    nu))

  sumlty = my_tofySICHEL2(
    as.double(x),
    as.double(mu),
    as.double(sigma),
    as.double(nu),
    as.double(lbes),
    as.double(cvec),
    ans = double(length(x)),
    as.integer(length(x)),
    as.integer(max(x) + 1)

  )


  sumlty <- as.double(
    .C(
      "tofySICHEL2",
      as.double(x),
      as.double(mu),
      as.double(sigma),
      as.double(nu),
      as.double(lbes),
      as.double(cvec),
      ans = double(length(x)),
      as.integer(length(x)),
      as.integer(max(x) + 1),
      PACKAGE = "gamlss.dist"
    )$ans
    )

  logfy <- -lgamma(x + 1) - nu * log(sigma * alpha) + sumlty +
    log(besselK(alpha, nu)) - log(besselK((1/sigma), nu))
  if (log == FALSE)
    fy <- exp(logfy)
  else fy <- logfy
  if (length(sigma) > 1)
    fy <- ifelse((sigma > 10000) & (nu > 0), dNBI(x, mu = mu,
                                                  sigma = 1/nu, log = log), fy)
  else fy <- if ((sigma > 10000) & (nu > 0))
    dNBI(x, mu = mu, sigma = 1/nu, log = log)
  else fy
  fy
}

my_dSICHEL(700, 100, 0.1, -30)
