/*  Variable naming:
  obs       = observed
  cen       = (right) censored
  N         = number of samples
  tau       = scale parameter
  J         = shared frailty parameter
  
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
    }
  
  data {
    int<lower=0> Nobs;
    int<lower=0> Ncen;
    vector[Nobs] yobs;
    vector[Ncen] ycen;
    int<lower=0> M;
    matrix[Nobs, M] Zobs;
    matrix[Ncen, M] Zcen;
  }
  
  transformed data {
    real<lower=0> tau_al;
    real<lower=0> tau_mu;
    int<lower=0> N;
    tau_al = 10.0;
    tau_mu = 10.0;
    N = Nobs + Ncen;

  }
  
  parameters {
    real alpha_raw;
    real mu;
    
    real<lower=0> tau_s_b_raw;
    vector<lower=0>[M] tau_b_raw;
    vector[M] beta_b_raw;
  }
  
  transformed parameters {
    real<lower=0> alpha;
    vector[M] beta_b;

    beta_b = b_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;
    alpha = exp(tau_al * alpha_raw);

  }
  
  model {
    //prior   
    beta_b_raw ~ normal(0.0, 1.0);
    alpha_raw ~ normal(0.0, 1.0);
    mu ~ normal(0 , tau_mu);
    
    //model
    yobs ~ weibull(alpha, exp(-( mu +  Zobs * beta_b )/alpha));
    target += weibull_lccdf(ycen | alpha, exp(-( mu +  Zcen * beta_b)/alpha));    
  }
  
  generated quantities {
    vector[N] log_lik;
    
    for (n in 1:Nobs){
      log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu + Zobs[n,] * beta_b)/alpha));
    }
    for (n in 1:Ncen){
      log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu + Zcen[n,] * beta_b)/alpha));
    }
    
  }
  
  
  