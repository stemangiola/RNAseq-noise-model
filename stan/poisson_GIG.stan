functions{
  real log_besselk_frac(real v, real z);

  real sichel_lpmf(int y, real gamma, real alpha, real epsilon){
    // https://www.researchgate.net/publication/247109928_Parameter_Estimation_for_the_Sichel_Distribution_and_Its_Multivariate_Extension


    real w = hypot(epsilon, alpha) - epsilon;

    return
      ( (gamma * (log(w) - log(alpha))) - log_besselk_frac(gamma, w) ) +
      (y * (log(epsilon) + log(w) - log(alpha))) - lgamma(y+1) +
      log_besselk_frac(y + gamma, alpha);

  }

}
data{
  int<lower=0> N;
  int y[N];

}
parameters{
  real<lower=-0.5, upper=0.5> gamma;
  real<lower=0, upper=0.5> alpha;
  real<lower=0, upper=0.5> epsilon;
}
model {
gamma ~ normal(0,1);
alpha ~ normal(0,1);
epsilon ~ normal(0,1);

// //nu ~ student_t(3, 0, 1);
// sigma ~ student_t(3, 0, 10);
print(gamma, " " , alpha, " " , epsilon);
for(i in 1:N) print(sichel_lpmf(y[i] | gamma, alpha , epsilon));
//for(i in 1:N) y[i] ~ sichel(gamma, alpha , epsilon);


//  print(sichel_lpmf(x | -30.0, 10.0 , 0.1));
//  print(sichel_lpmf(1000 | -2.0, 2.0 , 0.257));



}

