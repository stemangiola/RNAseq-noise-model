library(rstan)

besselk_test_stan <-
  '
functions {
  real log_besselk_frac(real v, real z);
  real log_besselk_frac_wrapper(real v, real z) {
    return log_besselk_frac(v,z);
  }
}

'

besselk_test_compiled <- stanc(model_code = besselk_test_stan, model_name = "besselk_test",
                                    allow_undefined = TRUE)

#cat(besselk_test_compiled$cppcode)

expose_stan_functions(besselk_test_compiled,
                      includes = paste0('\n#include "',here::here("dev","besselk.hpp"),'"\n'),
                      rebuild = TRUE)

#tt <- stan_model(model_code = besselk_test_stan, allow_undefined = TRUE,
#                 includes = paste0('\n#include "',here::here("dev","besselk.hpp"),'"\n'))

log_besselk_frac_wrapper(1.81,2.66)
