/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
 J1        = indicator to each subgroup

author :Carlos Traynor.
Grateful to Tomi Peltola and Jacqueline Buros for guidance
*/
data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  vector[Nobs] yobs;
  vector[Ncen] ycen;

}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  
  N = Nobs + Ncen;

  tau_al = 10.0;
  tau_mu = 10.0;
}

parameters {
  real alpha_raw;
  real mu;
  
}

transformed parameters {
  real<lower=0> alpha;

  alpha = exp(tau_al * alpha_raw);

}

model {
  yobs ~ weibull(alpha,  exp(-(mu)/alpha));
  target += weibull_lccdf(ycen | alpha,  exp(-(mu)/alpha));

  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);

}
generated quantities {
  vector[N] log_lik;
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-(mu)/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-(mu)/alpha));
  }
}
