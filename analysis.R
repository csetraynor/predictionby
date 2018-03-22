#Analysis
library(loo)
library(readr)
library(rstan)
library(rstanarm)
library(bayesplot)
library(rethinking)
theme_set(theme_bw())
stannull <- read_rds("stanfit/stannul.rds")
log_liknull <- loo::extract_log_lik(stannull, parameter_name = "log_lik")
colnames(log_liknull) <- md$intclust
loonull <- loo::loo(tbl_df(log_liknull))
#rm(list = c('stannull', 'log_liknull'))

stanclinical <- read_rds("stanfit/clin.rds")
log_likclinical <- loo::extract_log_lik(stanclinical, parameter_name = "log_lik")
looclinical <- loo::loo(log_likclinical)
#rm(list = c('stanclinical', 'log_likclinical'))
elpdclin = compare(looclinical, loonull)

stangenomic <- read_rds("stanfit/genomic.rds")
log_likgenomic <- loo::extract_log_lik(stangenomic, parameter_name = "log_lik")
loogenomic <- loo::loo(log_likgenomic)
#rm(list = c('stangenomic', 'log_likgenomic'))
elpdgene = compare(loogenomic, loonull)

colnames(log_likclinical) <- md$intclust
colnames(log_likgenomic) <- md$intclust
#loo estimate function
looestimate = function(x, log_lik){
  out <- loo::loo(log_lik[, colnames(log_lik) == x])
  return (out)
}
looestimate(y[[1]], log_likclinical)

#Int Clust
y <- as.list(c(1:3, "4ER+", "4ER-", 5:10))
# lolik <- list(c(log_liknull, log_likclinical, log_likgenomic))
#

## Create a dataframe with differences in leave one out cross validation respect to the null model stratified by clusters
null=lapply(y, function(x) looestimate(x, log_liknull))
clin=lapply(y, function(x) looestimate(x, log_likclinical))
gene=lapply(y, function(x) looestimate(x, log_likgenomic))
elpdclin = mapply(function(x, y) compare(x, y), x=clin, y=null)
elpdgene = mapply(function(x, y) compare(x, y), x=gene, y=null)


d = data.frame(Model = c(rep("C",11), rep("G",11)),
               IntClust = c(rep(y %>% unlist, 2)),
               elpd_diff = c(elpdclin[1,],elpdgene[1,]),
               elpd_diff_upper = c(elpdclin[1,] + 1.96*elpdclin[2,], elpdgene[1,] + 1.96*elpdgene[2,]),
               elpd_diff_lower = c(elpdclin[1,] - 1.96*elpdclin[2,], elpdgene[1,] - 1.96*elpdgene[2,]))
d$Model <- factor(d$Model, levels = c("Clinical", "Genomic Multilevel"))
d$IntClust = factor(d$IntClust, levels = c("1", "2", "3", "4ER+", "4ER-", "5",
                                             "6", "7", "8", "9", "10"))
ggplot2::ggplot(data = d , aes(x = Model, y = elpd_diff, ymin = elpd_diff_lower, ymax = elpd_diff_upper)) + geom_errorbar(width=0.2, size= .5, color="blue")  + geom_point( mapping=aes(x=Model, y=elpd_diff), size=2, shape=21, fill="white") + ylim(c( -40, 10)) + hline_0()+labs(title = "Leave One Out cross-validation error", subtitle = "Comparision with Clinical model by IntClust", y = "Leave One Out cross-validation error")+facet_wrap(~IntClust)

posterior <- read_rds("stanfit/genomic.rds")
posterior <- as.array(posterior, pars = "beta_p")

n_gene <- rownames(exprs(gd))
iclust <- as.list(c(1:3, "4ER+", "4ER-", 5:10))
names = lapply(iclust, function(x) rep(x, length = n_distinct(n_gene)))
dimnames(posterior)$parameters <- names %>% unlist

#Select which genes are have a larger absolute value posterior median 
select_by_iclust = function(x, level = .9){
  lower_p <- 0 + ((1 - level)/2)
  upper_p <- 1 - ((1 - level)/2)
  
  ppiclust <-  posterior[,,dimnames(posterior)$parameters == x]
  select <-  c(apply(ppiclust, 3, median) < quantile(apply(ppiclust, 3, median), probs = lower_p)|
               apply(ppiclust, 3, median) > quantile(apply(ppiclust, 3, median), probs = upper_p))
  dimnames(ppiclust)$parameters <- n_gene
  out <-  dimnames(ppiclust[,,select %>% unlist])$parameters
  return(out)
}

select = lapply(iclust, function(x) select_by_iclust(x, level= .9))
test_gene = unique(select %>% unlist)
test_gd = gd[featureNames(gd) %in% test_gene,]
saveRDS(test_gd, file = "gd_test.rds")
test_gd <- read_rds("gd_test.rds")
#central posterior uncertainty intervals


#Plot the marginal posterior distributions
posterior <- read_rds("stanfit/genomic.rds")
posterior <- as.array(posterior, pars = "beta_p")
x <- rownames(gd)
y <- as.list(c(1:3, "4ER+", "4ER-", 5:10))
names <- mapply(function(x) mapply( function(y)
  paste(rownames(gd), "IntClust", y , sep = "_"), y = y) , x = x)
dimnames(posterior)$parameters <- names[,1]
rm(names)
mcmc_hist(posterior, pars = dimnames(posterior)$parameters[grepl(paste0("(", paste(gene_select[[10]], collapse = "|"), ").*IntClust_10$"), dimnames(posterior)$parameters)]) +
  ggplot2::labs(
    title = "Posterior distributions",
    subtitle = "with medians, 50% and 90% intervals",
    y = "Shared Frailty parameters"
  )+ vline_0()


#uncertainty intervals as shaded areas
mcmc_areas(
  posterior, 
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)


