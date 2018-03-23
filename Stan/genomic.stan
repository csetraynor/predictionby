/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
 G         = gmups#
 mu        = shared frailty var
 

author Carlos Traynor
*/
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }
    return res;
  }

  vector b_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
    return r_global * sqrt_vec(r_local);
  }
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu_local, real nu_global, real scale_global) {
    //half-t preior for lambdas
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
    //half-t prior for tau
    r1_global ~ normal(0.0, scale_global);
    r2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> G; //intClust
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int Jobs[Nobs]; 
  int Jcen[Ncen]; 
  int<lower=0> M;
  matrix[Nobs, M] Zobs;
  matrix[Ncen, M] Zcen;
  int<lower=0> P;
  matrix[Nobs, P] Pobs;
  matrix[Ncen, P] Pcen;
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  int<lower=0> nu_global;
  real<lower=0> nu_local;
  real<lower=0> scale_global;
  tau_al = 10.0;
  tau_mu = 10.0;
  N = Nobs + Ncen;
  nu_global = 100;
  nu_local = 1;
  scale_global = 1;
}

parameters {
  real alpha_raw;
    
  real<lower=0> tau_s_b_raw;
  vector<lower=0>[M] tau_b_raw;
  vector[M] beta_b_raw;
  
  real<lower=0> tau_s1_p_raw;
  real<lower=0> tau_s2_p_raw;
  vector<lower=0>[P] tau_1_p_raw;
  vector<lower=0>[P] tau_2_p_raw;
  matrix[P, G] beta_p_raw;

  real mu;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;
  matrix[P, G] beta_p;
  vector[Nobs] lpobs;
  vector[Ncen] lpcen;
  
  alpha = exp(tau_al * alpha_raw);
  beta_b = b_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;
  for(g in 1:G){
     beta_p[, g] = hs_prior_lp(tau_s1_p_raw, tau_s2_p_raw, tau_1_p_raw, tau_2_p_raw, nu_local, nu_global, scale_global) .* beta_p_raw[, g]; 
  }
  for (n in 1:Nobs){
    lpobs[n] = mu + Pobs[n,] * beta_p[, Jobs[n]] + Zobs[n,] * beta_b;
  }
  for (n in 1:Ncen){
    lpcen[n] = mu + Pcen[n,] * beta_p[, Jcen[n]] + Zcen[n,] * beta_b;
  }
}
model {
  //priors
  alpha_raw ~ normal(0.0, 1.0);
  
  beta_b_raw ~ normal(0.0, 1.0);
  
  for(g in 1:G){
    beta_p_raw[,g] ~ normal(0.0, 1.0); 
  }
  mu ~ normal(0 , tau_mu);
  
  //model
  yobs ~ weibull(alpha, exp(-(lpobs)/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-(lpcen)/alpha));
}
generated quantities {
  vector[N] yhat_uncens;
  vector[N] log_lik;
  
  //log likelihood
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-(lpobs[n])/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-(lpcen[n])/alpha));
  }
  //yhat uncens calculation
  for (n in 1:Nobs){
    yhat_uncens[n] = weibull_rng(alpha, exp(-(lpobs[n])/alpha));
  }
  for (n in 1:Ncen){
    yhat_uncens[Nobs + n] = weibull_rng(alpha, exp(-(lpcen[n])/alpha));
  }
  
}

