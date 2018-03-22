# !diagnostics off
rm(list = ls())
library(readr)
library(dplyr)
library(caret)
library(Biobase)

##--- Gen Data Null
gen_stan_data0 <- function(data) {

  ind = (data$status == 1)
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind]
  )
}
gen_inits0 <- function() {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1)
    )
}
 # into_data <- gen_stan_data(md)
 # glimpse(into_data)
 # rstan::stan_rdump(ls(into_data), file = "checking.data.R",
 #                   envir = list2env(into_data))



##--- Gen Stan Clinical Model ---#
gen_stan_data1 <- function(data, formula = as.formula(~1)) {
  #Covariates (clinical) Matrix
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- as.matrix(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }
  
  ind = data$status == 1
  Zobs <- Z[ind,]
  Zcen <- Z[!ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[!ind,]), M))
  )
}
gen_stan_data1(md , formula = "~ tumor_stage + size + grade ") %>%  glimpse
# glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits1 <- function(M) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      beta_b_raw = rnorm(M)
    )
}

# inits <- gen_inits(J = 6)
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))

##--- Gen Stan Multilevel ClinicoGenomic Model ---#
gen_stan_data2 <- function(data, formula = as.formula(~1), group = "intclust", newvar = eset) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  #Covariates (clinical) Matrix
  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- as.matrix(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }

  # Subgroup indicator
  Jgroup <- data%>% select(group) %>% unlist %>% as.factor %>% as.integer() 
  
  Peta <- t(exprs(newvar))
  I_P = ncol(Peta)
  ind = data$status == 1
  Zobs <- Z[ind,]
  Zcen <- Z[!ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[!ind,]), M)),
    G = n_distinct(Jgroup),
    Jobs = Jgroup[ind],
    Jcen = Jgroup[!ind],
    P = I_P,
    Pobs = Peta[ind,],
    Pcen = Peta[!ind,]
  )
}


gen_stan_data2(md, formula = "~ tumor_stage + size + grade", newvar = gd) %>%  glimpse

#   glimpse

# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits2 <- function(J, M, P) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      tau_s1_p_raw = 0.1*abs(rnorm(1)),
      tau_s2_p_raw = 0.1*abs(rnorm(1)),
      tau_1_p_raw = abs(rnorm(P)),
      tau_2_p_raw = abs(rnorm(P)),
      
      beta_b_raw = rnorm(M),
      beta_p_raw = array(rnorm(P*J), dim = c(P, J)),
      
      mu = rnorm(1),
      
      peta = rnorm(P),
      peta_sigma = abs(rnorm(P))
      
    )
}
# init
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))
gen_inits3 <- function(J, M, P) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      tau_s_g_raw = 0.1*abs(rnorm(1)),
      tau_g_raw = abs(rnorm(P)),
      
      beta_b_raw = rnorm(M),
      beta_p_raw = array(rnorm(P*J), dim = c(P, J)),
      
      mu = rnorm(1),
      
      peta = rnorm(P),
      peta_sigma = abs(rnorm(P))
      
    )
}


save(list = ls(), file = "Gen_data_fun.Rdata")
rm(list = ls())
