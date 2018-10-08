#logsumexp + softmax from https://gist.github.com/aufrank/83572

logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

#' Simulate resampling count data - given counts of various items, each item is given
#' equal probability to be sampled.
#' @param draw A vector of item counts
#' @param k Number of items to be selected
resample_draw = function(draw, k) {
  if(k >= sum(draw)) {
    stop("The number reads to sample (k) has to be smaller than the sum of draw")
  }
  G = length(draw)
  positions = c(0, cumsum(draw))
  read_set = array(-1, sum(draw))
  for(g in 1:G) {
    read_set[(positions[g] + 1):positions[g + 1]] = g
  }

  sampling_result = sample(read_set, k)

  counts = array(-1, G)
  for(g in 1:G) {
    counts[g] = sum(sampling_result == g)
  }

  counts
}

#' Generate a vector of size N that sums to zero and marginal
#' distributions for all elements are N(0,1)
#' Original code due to andre.pfeuffer
#' https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884/31
rnorm_sum_to_zero = function(N) {
  Q_r = numeric(2 * N)
  scale = 1/sqrt(1-1/N)
  for(i in 1:N) {
    Q_r[i] = -sqrt((N-i)/(N-i + 1)) * scale
    Q_r[i+N] = 1/sqrt((N - i) * (N - i + 1)) * scale
  }

  x_raw = rnorm(N - 1, 0, 1)
  x = numeric(N)
  x_aux = 0

  for(i in 1:(N-1)) {
    x[i] = x_aux + x_raw[i] * Q_r[i]
    x_aux = x_aux + x_raw[i] * Q_r[i + N]
  }
  x[N] = x_aux
  x
}
