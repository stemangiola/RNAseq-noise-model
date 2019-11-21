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

waring_lpdf <- function(y, mu, k, rho) {
  a <- mu * (rho - 1) / k
  lgamma(a + rho) +
    lgamma(k + rho) -
    lgamma(rho) -
    lgamma(a + k + rho) +
    lgamma(a + y) -
    lgamma(a) +
    lgamma(k + y) -
    lgamma(k) -
    lgamma(a + k + rho + y) +
    lgamma(a + k + rho) -
    lgamma(y + 1)
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


