functions{
  real log_besselk_frac(real v, real z);
}
parameters{
  real alpha;
}
model {
  alpha ~ normal(0,1);
  print(log_besselk_frac(0.1, 100) );
  print(log_besselk_frac(10, 100) );
  print(log_besselk_frac(30, 100) );
  print(log_besselk_frac(30, 1000) );

}

