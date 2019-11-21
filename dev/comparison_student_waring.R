library(tidyverse)

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 12),
		aspect.ratio = 1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.title.y  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		))
	)

# STUDENT T

main_data <-

	tibble(mu = 50) %>%
	crossing(tibble(y = seq(from = -250, to = 250, by = 1))) %>%
	crossing(tibble(phi = c( 5)))%>%
	crossing(tibble(tau = c( seq(1, 10, 1), 100, 1000))) %>%
	mutate(		ld = LaplacesDemon::dst(y, mu = mu, sigma = phi, nu = tau, log = TRUE) ) %>%
	bind_rows(
		tibble(mu = 500) %>%
			crossing(tibble(y = seq(from = -2000, to = 2000, by = 10))) %>%
			crossing(tibble(phi = c(50)))%>%
			crossing(tibble(tau = c( seq(1, 10, 1), 100, 1000)))%>%
			mutate(		ld = LaplacesDemon::dst(y, mu = mu, sigma = phi, nu = tau, log = TRUE) )
	) %>%
	bind_rows(
		tibble(mu = 5000) %>%
			crossing(tibble(y = seq(from = -20000, to = 20000, by = 100))) %>%
			crossing(tibble(phi = c( 500))) %>%
			crossing(tibble(tau = c(seq(1, 10, 1), 100, 1000))) %>%
			mutate(		ld = LaplacesDemon::dst(y, mu = mu, sigma = phi, nu = tau, log = TRUE) )
	) %>%
	bind_rows(
		tibble(mu = 10000) %>%
			crossing(tibble(y = seq(from = -50000, to = 50000, by = 100))) %>%
			crossing(tibble(phi = c( 1000)))%>%
			crossing(tibble(tau = c(seq(1, 10, 1), 100, 1000))) %>%
			mutate(		ld = LaplacesDemon::dst(y, mu = mu, sigma = phi, nu = tau, log = TRUE) )
	)

normal_data <-
	main_data %>%
	distinct(mu ,    y ,  phi) %>%
	mutate(
		nb_dens = dnorm(y, mean = mu, sd = phi, log = TRUE)
	)

# Relationship tau mu LOWER
main_data %>%
	distinct( mu,  phi,   tau ) %>%
	mutate(q2.5 = LaplacesDemon::qst(0.05, mu, phi, tau)) %>%
	ggplot(aes(x = tau, y = q2.5/mu, color = interaction(mu, phi))) + geom_jitter() + scale_x_log10() + my_theme

# Relationship tau mu UPPER
main_data %>%
  distinct( mu,  phi,   tau ) %>%
  mutate(q95 = LaplacesDemon::qst(0.95, mu, phi, tau)) %>%
  ggplot(aes(x = tau, y = q95/mu, color = interaction(mu, phi))) + geom_jitter() + scale_x_log10() + my_theme

(main_data %>%
		ggplot(aes(x = y, y = ld %>% exp, color = tau)) +
		geom_line(aes(group = tau)) +
		geom_line(data = normal_data, aes(y = nb_dens %>% exp), color="red") +
		facet_wrap(~mu+phi, ncol = 2, labeller = "label_both", scales = "free") +
		my_theme) %>% plotly::ggplotly()


# tau pathological behaviour

main_data <-
	tibble(mu = 5) %>%
	crossing(tibble(y = seq(from = -50, to = 50, by = 1))) %>%
	crossing(tibble(phi = 1))%>%
	crossing(tibble(tau = c( seq(1, 10, 1), seq(20, 100, 10), seq(200, 2000, 200)))) %>%
	mutate(		ld = LaplacesDemon::dst(y, mu = mu, sigma = phi, nu = tau, log = TRUE) )

nb_data <-
	main_data %>%
	distinct(mu ,    y ,  phi) %>%
	mutate(
		nb_dens = dnorm(y, mean = mu, sd = phi, log = TRUE)
	)

(main_data %>%
		ggplot(aes(x = y, y = ld , color = tau)) +
		geom_line(aes(group = tau)) +
		geom_line(data = nb_data, aes(y = nb_dens ), color="black") +
		facet_wrap(~mu+phi,   scales = "free") +
		my_theme) %>% plotly::ggplotly()

# WARING

gwaring_lpmf <- function(y, mu, k, rho) {
  a <- mu * (rho - 1) / k
  lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) -
    lgamma(a + k + rho) + lgamma(a + y) - lgamma(a) +
    lgamma(k + y) - lgamma(k) - lgamma(a + k + rho + y) +
    lgamma(a + k + rho) - lgamma(y + 1)
}


waring_lpdf <- function(y, mu, k, rho) {
	a <- mu * (rho - 1) / k
	lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) -
	  lgamma(a + k + rho) + lgamma(a + y) - lgamma(a) +
	  lgamma(k + y) - lgamma(k) - lgamma(a + k + rho + y) +
	  lgamma(a + k + rho) - lgamma(y + 1)
}

rwaring <- function(n, mu, phi, tau) {
  k <-  (tau + 1)*phi
  rho <-  (phi + 1) / tau + phi + 2
  a <-  mu * (rho - 1) / k

  v_logit <- rbeta(n, k, rho)

  v <-  v_logit / (1 - v_logit)
  neg_bin_mu = a * v
  rnbinom(n, mu = neg_bin_mu, size = a)
}


# Fixing the decomposition of variance
# https://www.wolframalpha.com/input/?i=solve+for+k%2C+rho+%3A+tau+%3D+%28k+%2B+1%29%2F%28rho+-+2%29%3B+phi+%3D+k+*+%28rho+-+2%29+%2F+%28k+%2B+rho++-1%29

# tau not correctly scaled

main_data <-

	tibble(mu = 50) %>%
	crossing(tibble(y = seq(from = 0, to = 1000, by = 1))) %>%
	crossing(tibble(phi = 50))%>%
	crossing(tibble(tau = c(0.001, 0.01, seq(0.1, 1, 0.1), seq(1, 10, 1), 100, 1000))) %>%
	mutate(
		k = (tau + 1)*phi,
		rho = (phi + 1) / tau + phi + 2,
		ld = waring_lpdf(y, mu = mu, k = k , rho = rho)) %>%
	bind_rows(
		tibble(mu = 500) %>%
			crossing(tibble(y = seq(from = 0, to = 5000, by = 10))) %>%
			crossing(tibble(phi = 50))%>%
			crossing(tibble(tau = c(0.001, 0.01, seq(0.1, 1, 0.1), seq(1, 10, 1), 100, 1000))) %>%
			mutate(
				k = (tau + 1)*phi,
				rho = (phi + 1) / tau + phi + 2,
				ld = waring_lpdf(y, mu = mu, k = k , rho = rho))
	) %>%
	bind_rows(
		tibble(mu = 5000) %>%
			crossing(tibble(y = seq(from = 0, to = 50000, by = 100))) %>%
			crossing(tibble(phi = 50))%>%
			crossing(tibble(tau = c(0.001, 0.01, seq(0.1, 1, 0.1), seq(1, 10, 1), 100, 1000))) %>%
			mutate(
				k = (tau + 1)*phi,
				rho = (phi + 1) / tau + phi + 2,
				ld = waring_lpdf(y, mu = mu, k = k , rho = rho))
	) %>%
	bind_rows(
		tibble(mu = 10000) %>%
			crossing(tibble(y = seq(from = 0, to = 50000, by = 100))) %>%
			crossing(tibble(phi = 50))%>%
			crossing(tibble(tau = c(0.001, 0.01, seq(0.1, 1, 0.1), seq(1, 10, 1), 100, 1000))) %>%
			mutate(
				k = (tau + 1)*phi,
				rho = (phi + 1) / tau + phi + 2,
				ld = waring_lpdf(y, mu = mu, k = k , rho = rho))
	)

nb_data <-
	main_data %>%
	distinct(mu ,    y ,  phi) %>%
	mutate(
		nb_dens = dnbinom(y, mu = mu, size = phi, log = TRUE)
	)

# Relationship tau mu LOWER
main_data %>%
	distinct( mu,  phi, tau ) %>%
  rowwise() %>%
	mutate(q5 = rwaring(10000, mu, phi, tau) %>% quantile(c(0.05))) %>%
  ungroup %>%
	ggplot(aes(x = 1/tau, y = q5/mu, color = factor(mu))) + geom_jitter() + scale_x_log10() + my_theme

# Relationship tau mu UPPER
main_data %>%
  distinct( mu,  phi, tau ) %>%
  rowwise() %>%
  mutate(q95 = rwaring(10000, mu, phi, tau) %>% quantile(c(0.95))) %>%
  ungroup %>%
  ggplot(aes(x = 1/tau, y = q95/mu, color = factor(mu))) + geom_jitter() + scale_x_log10() + my_theme


(main_data %>%
		ggplot(aes(x = y, y = ld %>% exp, color = tau)) +
		geom_line(aes(group = tau)) +
		geom_line(data = nb_data, aes(y = nb_dens %>% exp), color="red") +
		facet_wrap(~mu+phi, ncol = 2, labeller = "label_both", scales = "free") +
		my_theme) %>% plotly::ggplotly()


# tau pathological behaviour

main_data <-
	tibble(mu = 5000) %>%
	crossing(tibble(y = seq(from = 0, to = 50000, by = 100))) %>%
	crossing(tibble(phi = c(50, 500, 5000)))%>%
	crossing(tibble(tau = c(5, 50, 500))) %>%
	mutate(
		k = (tau + 1)*phi,
		rho = (phi + 1) / tau + phi + 2,
		ld = waring_lpdf(y, mu = mu, k = k , rho = rho))

nb_data <-
	main_data %>%
	distinct(mu ,    y ,  phi) %>%
	mutate(
		nb_dens = dnbinom(y, mu = mu, size = phi, log = TRUE)
	)

(main_data %>%
		ggplot(aes(x = y, y = ld , color = factor(tau))) +
		geom_line(aes(group = tau)) +
		geom_line(data = nb_data, aes(y = nb_dens ), color="black") +
		facet_wrap(~mu+phi,   scales = "free") +
		my_theme) %>% plotly::ggplotly()


# WARING fixed
my_tau = c(0.001, 0.01, seq(0.1, 1, 0.1), seq(1, 10, 1))/100

waring_fixed_lpmf <- function(y, mu, phi, tau) {

  tau = tau * mu

  k = (tau + 1)*phi
  rho = (phi + 1) / tau + phi + 2
  waring_lpdf(y, mu = mu, k = k , rho = rho)

}

rwaring_fixed <- function(n, mu, phi, tau) {

  tau = tau * mu

  k <-  (tau + 1)*phi
  rho <-  (phi + 1) / tau + phi + 2
  a <-  mu * (rho - 1) / k

  v_logit <- rbeta(n, k, rho)

  v <-  v_logit / (1 - v_logit)
  neg_bin_mu = a * v
  rnbinom(n, mu = neg_bin_mu, size = a)
}

main_data <-

  tibble(mu = 50) %>%
  crossing(tibble(y = seq(from = 0, to = 1000, by = 1))) %>%
  crossing(tibble(phi = 50))%>%
  crossing(tibble(tau = my_tau)) %>%
  mutate(  ld = waring_fixed_lpmf(y, mu = mu, phi = phi, tau = tau)) %>%
  bind_rows(
    tibble(mu = 500) %>%
      crossing(tibble(y = seq(from = 0, to = 5000, by = 10))) %>%
      crossing(tibble(phi = 50))%>%
      crossing(tibble(tau = my_tau)) %>%
      mutate(  ld = waring_fixed_lpmf(y, mu = mu, phi = phi, tau = tau))
  ) %>%
  bind_rows(
    tibble(mu = 5000) %>%
      crossing(tibble(y = seq(from = 0, to = 50000, by = 100))) %>%
      crossing(tibble(phi = 50))%>%
      crossing(tibble(tau = my_tau)) %>%
      mutate(  ld = waring_fixed_lpmf(y, mu = mu, phi = phi, tau = tau))
  ) %>%
  bind_rows(
    tibble(mu = 10000) %>%
      crossing(tibble(y = seq(from = 0, to = 50000, by = 100))) %>%
      crossing(tibble(phi = 50))%>%
      crossing(tibble(tau = my_tau)) %>%
      mutate(  ld = waring_fixed_lpmf(y, mu = mu, phi = phi, tau = tau))
  )


nb_data <-
  main_data %>%
  distinct(mu ,    y ,  phi) %>%
  mutate(
    nb_dens = dnbinom(y, mu = mu, size = phi, log = TRUE)
  )

# Relationship tau mu LOWER
main_data %>%
  distinct( mu,  phi, tau ) %>%
  rowwise() %>%
  mutate(q5 = rwaring_fixed(10000, mu, phi, tau) %>% quantile(c(0.05))) %>%
  ungroup %>%
  ggplot(aes(x = 1/tau, y = q5/mu, color = factor(mu))) + geom_jitter() + scale_x_log10() + my_theme

# Relationship tau mu UPPER
main_data %>%
  distinct( mu,  phi, tau ) %>%
  rowwise() %>%
  mutate(q95 = rwaring_fixed(10000, mu, phi, tau) %>% quantile(c(0.95))) %>%
  ungroup %>%
  ggplot(aes(x = 1/tau, y = q95/mu, color = factor(mu))) + geom_jitter() + scale_x_log10() + my_theme


# Relationship tau mu
(main_data %>%
    ggplot(aes(x = y, y = ld %>% exp, color = tau)) +
    geom_line(aes(group = tau)) +
    geom_line(data = nb_data, aes(y = nb_dens %>% exp), color="red") +
    facet_wrap(~mu+phi, ncol = 2, labeller = "label_both", scales = "free") +
    my_theme) %>% plotly::ggplotly()




# Avoiding pathological tau lower tail

# wrong shape



tau = 0.03
mu = 10000
phi = 50
tau = tau * mu
k <-  (tau + 1)*phi
rho <-  (phi + 1) / tau + phi + 2
a <-  mu * (rho - 1) / k
rbeta(10000, k, rho) %>% density %>% plot

tibble(
w = rwaring_fixed(100000, mu, phi, tau) ,
n = rnbinom(100000, mu = mu, size = phi)
) %>%
  gather(dist, count) %>%
  ggplot(aes(count, color=dist)) + geom_density() + scale_y_log10()

main_data %>%
  filter(mu == !!mu & phi  == !!phi & tau == !!tau) %>%
  ggplot(aes(x = y, y = ld %>% exp, color = tau)) +
  geom_line(aes(group = tau)) +
  geom_line(
    data = nb_data %>% filter(mu == !!mu & phi  == !!phi & tau == !!tau),
    aes(y = nb_dens %>% exp), color="red") +
  facet_wrap(~mu+phi, ncol = 2, labeller = "label_both", scales = "free")
  my_theme

tau = 0.02
mu = 10000
phi = 50
tau = tau * mu
k <-  (tau + 1)*phi
rho <-  (phi + 1) / tau + phi + 2
a <-  mu * (rho - 1) / k
rbeta(10000, k, rho) %>% density %>% plot

tau = 0.01
mu = 10000
phi = 50
tau = tau * mu
k <-  (tau + 1)*phi
rho <-  (phi + 1) / tau + phi + 2
a <-  mu * (rho - 1) / k
rbeta(10000, k, rho) %>% density %>% plot

list(
  tibble(
    w = rwaring_fixed(100000, mu, phi, tau) ,
    n = rnbinom(100000, mu = mu, size = phi)
  ) %>%
    gather(dist, count) %>%
    ggplot(aes(count, color=dist)) + geom_density() + scale_y_log10() + scale_color_manual(values = c("red", "blue")),

  main_data %>%
    filter(mu == !!mu & phi  == !!phi & tau == !!tau) %>%
    ggplot(aes(x = y, y = ld %>% exp, color = tau)) +
    geom_line(aes(group = tau)) +
    geom_line(
      data = nb_data %>% filter(mu == !!mu & phi  == !!phi & tau == !!tau),
      aes(y = nb_dens %>% exp), color="red") + scale_y_log10() +
    facet_wrap(~mu+phi, ncol = 2, labeller = "label_both", scales = "free")
) %>% cowplot::plot_grid(plotlist = .)



tau = 0.001
mu = 10000
phi = 50
tau = tau * mu
k <-  (tau + 1)*phi
rho <-  (phi + 1) / tau + phi + 2
a <-  mu * (rho - 1) / k
rbeta(10000, k, rho) %>% density %>% plot

tau = 1e-5
mu = 10000
phi = 50
tau = tau * mu
k <-  (tau + 1)*phi
rho <-  (phi + 1) / tau + phi + 2
a <-  mu * (rho - 1) / k
rbeta(10000, k, rho) %>% density %>% plot
