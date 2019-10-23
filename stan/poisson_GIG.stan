functions{
  real log_besselk_frac(real v, real z);

  real[] tofySICHEL2(int[] y, real mu, real sigma, real nu, real lbes, real cvec){

  int maxyp1 = max(y) + 1;
  vector[maxyp1] tofY;
  int iy;
  real alpha;
  real sumT;
  real ans[size(y)];

  for (i in 1:size(y))  {
    iy = y[i]+1;
    tofY[1] = (mu/cvec)*pow((1.0+2.0*sigma*mu/cvec),(-0.5))*exp(lbes);
    alpha = sqrt(1.0 + 2.0*sigma*mu/cvec)/sigma;
    sumT = 0.0;

    for (j in 2:iy)  {
      tofY[j] = ( cvec*sigma*(2.0*((j-1.0)+nu)/mu) + (1.0/tofY[j-1])) * pow(mu/(sigma*alpha*cvec),2);
      sumT += log(tofY[j-1]);
    }
    ans[i] = sumT;
  }
  return ans;
}

  real sichel_lpmf(int[] y, real nu, real mu, real sigma){


      if (sigma > 10000 && nu > 0) return neg_binomial_2_log_lpmf(y | mu, 1.0/nu);
      else{

        real nu_1 = nu + 1.0;
        real one_over_sigma = 1.0/sigma;
        real log_bes_nu_one_over_sigma = log_besselk_frac(nu, one_over_sigma);
        real cvec = exp(log_besselk_frac(nu_1, one_over_sigma) - log_bes_nu_one_over_sigma);
        real alpha = sqrt(1.0 + 2.0 * sigma * (mu/cvec))/sigma;
        real log_bes_nu_alpha = log_besselk_frac(nu, alpha);
        real lbes = log_besselk_frac(nu_1, alpha) - log_bes_nu_alpha;
        real sumlty[size(y)] = tofySICHEL2(y, mu,sigma,nu,lbes,cvec);

        real lp[size(y)];

//print(sumlty);
       for(i in 1:size(y))
          lp[i] =  -lgamma(y[i] + 1) - nu * log(sigma * alpha) + sumlty[i] + log_bes_nu_alpha - log_bes_nu_one_over_sigma;

        return sum(lp);

      }
  }

}
data{
  int<lower=0> N;
  int y[N];

}
parameters{
  real<upper=0> nu;
  real<lower=0> mu;
  real<lower=0> sigma;
}
model {
  int x[1];
x[1] = 10;

//nu ~ student_t(3, 0, 1);
sigma ~ student_t(3, 0, 10);
y ~ sichel_lpmf(nu, mu , sigma);

  // print(sichel_lpmf(x | -30.0, 10.0 , 0.1));
  //   print(sichel_lpmf(x | 30.0, 10.0 , 0.1));

}
