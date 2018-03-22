##Run Stan
# !diagnostics off
library(dplyr)
library(readr)
library(survival)
library(rstan)
library(caret)
gd <- read_rds("Gen_Data.rds")
md <- read_rds("Med_Data_Clean.rds")
load("Gen_data_fun.Rdata")
# Run null model
stan_file_null <- "stan/null.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stannull <- rstan::stan(stan_file_null,
                        data = gen_stan_data0(md),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits0(),
                        control = list(adapt_delta = 0.95))
log_liknull <- loo::extract_log_lik(stannull, parameter_name = "log_lik")
loonull <- loo::loo(log_liknull)
print(loonull)
saveRDS(stannull, file = "stanfit/stannul.rds")
rm(list = c('stannull', 'log_liknull'))


#Run Clinical Model
stan_file_clin<- "stan/clinical.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanclin <- rstan::stan(stan_file_clin,
                        data = gen_stan_data1(md,
                        formula = "~ size + grade + tumor_stage"),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits1(M = 5),
control = list(adapt_delta = 0.95))
log_likclin <- loo::extract_log_lik(stanclin, parameter_name = "log_lik")
looclin <- loo::loo(log_likclin)
print(looclin)
compare(loonull, looclin) #preference for the second model!
saveRDS(stanclin, file = "stanfit/clin.rds")
rm(list = c('stanclin', 'log_likclin'))
stanclin=read_rds("stanfit/clin.rds")
#---------------------------------------
##Run Multilevel Genomic and Clinical Model for Selection
stan_file_genomic <- "stan/genomic.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stangenomic <- rstan::stan(stan_file_genomic,
                        data = gen_stan_data2(data = md,
                                              formula = "~ size + grade + tumor_stage",
                                              newvar = gd),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits2(J = 11, M = 5, P = 1219),
                        control = list(adapt_delta = 0.95))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

log_likgenomic <- loo::extract_log_lik(stangenomic, parameter_name = "log_lik")
loogenomic <- loo::loo(log_likgenomic)
print(loogenomic)
loo::compare(looclin, loogenomic)
saveRDS(stangenomic, file = "stanfit/genomic.rds")
rm(list = c('stangenomic', 'log_likgenomic'))
stangenomic=read_rds("stanfit/genomic.rds")





