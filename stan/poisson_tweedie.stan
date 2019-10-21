functions{

  vector zhuprobsLog(int n, real a, real b, real c){
    vector[n+1] res_log;
    vector[n] r_log ;
    if(a==0) res_log[1] = b* log(1-c);
    else res_log[1] = (b)*(((1.0-c)^a)-1.0)/(a);

    if(n!=0){
      real aux_log = log(b)+log(c);
      res_log[2] = aux_log+res_log[1];


      if(n > 1){
        r_log[1] = log(1-a)+log(c);
        for(i in 1:(n-1)) r_log[i+1] = log(c)+r_log[i]+log(i-1+a)-log(i+1);

        for(i in 1:(n-1)){
          res_log[i+1+1] = aux_log+res_log[i+1];
          for(j in 1:i) res_log[i+1+1] = log_sum_exp( res_log[i+1+1] , log(j)+r_log[i-j+1]+res_log[j+1]);
          res_log[i+1+1] = res_log[i+1+1] - log(i+1);
        }
      }
    }

    return res_log;
  }

  real poisson_tweedie_log_lpmf (int[] x, real log_mu, real D, real a) {

     real b = (exp(log_mu) * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a));
     real c = (D - 1)/(D - a);
     int xp1[size(x)];

     for(i in 1:size(x)) xp1[i] = x[i] + 1;
     return sum(zhuprobsLog(max(x), a, b, c)[xp1]);

   }

  real poisson_tweedie_raw(int[] x, real a, real b, real c);
  real poisson_tweedie_raw_da(int[] x, real a, real b, real c);

  real poisson_tweedie_raw_wrapper (int[] x, real log_mu, real D, real a) {
     real b = (exp(log_mu) * (1 - a)^(1 - a))/((D - 1) * (D - a)^(-a));
    real c = (D - 1)/(D - a);
    real xx = 1;
    return poisson_tweedie_raw(x, a, b, c);
  }

  real poisson_tweedie_abc_wrapper (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw(x, a, b, c);
  }

  real poisson_tweedie_abc_da_wrapper (int[] x, real a, real b, real c) {
    return poisson_tweedie_raw_da(x, a, b, c);
  }

}
data{
  int<lower=0> N;
  int y[N];

}
parameters{
  real mu;
  real<lower=0> sigma;
  real<lower=0, upper=1> nu;
}
model {
  mu ~ normal(0,2);
  sigma ~ student_t(3, 0, 1);
  nu ~ beta(1,5);
  y ~ poisson_tweedie_log(mu, sigma, nu);

}

