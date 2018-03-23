#Frailty Gen Data
##--- Gen Stan Multilevel ClinicoGenomic Model ---#
library(readr)
library(Biobase)
library(dplyr)

gd <- read_rds("Gen_Data.rds")
md <- read_rds("Med_Data_Clean.rds")
short_md <- md %>% dplyr::filter(intclust == "7")


# estimate the global parameter for shrinkage
p0 = 10  #choice of non zero parameters
D = nrow(exprs(gd))
n = nrow(short_md) 
library(MASS)
?fitdistr
fit=fitdistr(short_md %>% filter(status == 1) %>% dplyr::select(time) %>% unlist, densfun = "weibull")
fit$sd

scale_global = (p0 / (D - p0)) * (fit$sd[1] / (n)^(1/2)) 


gen_stan_data <- function(data, formula = as.formula(~1), group = "cohort", newvar = eset, scale_global = scale_global) {
  
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
  Jgroup <- data%>% dplyr::select(group) %>% unlist %>% as.factor %>% as.integer() 
  
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
    Pcen = Peta[!ind,],
    scale_global = scale_global
  )
}

gen_stan_data(short_md, formula = "~ tumor_stage + size + grade", newvar = gd, scale_global = scale_global ) %>%  glimpse

# prepare data for stan
into_data <- gen_stan_data(short_md, formula = "~ tumor_stage + size + grade", newvar = gd, scale_global = scale_global )
rstan::stan_rdump(ls(into_data), file = "checking.data.R",
                   envir = list2env(into_data))

gen_inits <- function(J, M, P) {
 # function()
    list(
      alpha_raw = 0.01*rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      tau_s1_p_raw = 0.1*abs(rnorm(1)),
      tau_s2_p_raw = 0.1*abs(rnorm(1)),
      tau_1_p_raw = abs(rnorm(P)),
      tau_2_p_raw = abs(rnorm(P)),
      
      beta_b_raw = rnorm(M),
      beta_p_raw = rnorm(P),
      
      mu = rnorm(J),
      tau_mu = abs(rnorm(1)),
      xi = rnorm(1)
      
      
    )
}
init = gen_inits(J = 5, M = 5, P = 2438)
rstan::stan_rdump(ls(init), file = "checking.init.R",
                                    envir = list2env(init))


