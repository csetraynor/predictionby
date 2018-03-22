#Data preprocess
library(readr)
library(VIM)
library(dplyr)
#load medical data
#md <- tbl_df(read_rds("C:/RFactory/bymetabric_files/rdsmetabric/medicaldfs.rds"))
md <- tbl_df(read_rds("medicaldfs.rds"))

#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    return(x)
    print("Error ", x, " not character")
  } else {
    ifelse(x == "" | x == "[Not Available]" | x == "[Discrepancy]" | x == "Indeterminate"| x == "Equivocal", NA, x)
  }
}
md <- md %>% dplyr::mutate_all(funs(convert_blank_to_na))
md %>% 
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
                 sortVars = TRUE, sortCombs = TRUE, plot = TRUE, only.miss = FALSE)

md$grade[md$lymph_nodes_positive == 0] <- "1"
md$grade[md$lymph_nodes_positive >=1 & md$lymph_nodes_positive <= 4 ] <- "2"
md$grade[md$lymph_nodes_positive > 4] <- "3"
md$lymph_nodes_positive <- NULL
#standardise continuos covariates and impute with median value
preProc <- caret::preProcess(md %>% select(size) %>% as.matrix(), method = c("center", "scale"))
md$size <- predict(preProc, md %>% select(size) %>% as.matrix())
hist(md$size)
md$size[is.na(md$size)] <- median(md$size, na.rm = T)
#Filter missing intclust
md <- md %>% filter(!is.na(intclust))
md %>% 
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
            sortVars = TRUE, sortCombs = TRUE, plot = TRUE, only.miss = FALSE)

#Imputation using Bayesian Logistic Regression
tmp <- as.factor(md$grade)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~ size +
                                    intclust  ,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$grade[is.na(md$grade)] <- tmp[is.na(md$grade)]
remove(tmp)
table(md$tumor_stage)
md$tumor_stage[md$tumor_stage == 0] <- 1
md$tumor_stage[md$tumor_stage == 4] <- 3

VIM::marginplot(md[c("time","tumor_stage")])
table(md$tumor_stage , md$intclust, useNA = "always")
tmp <- as.factor(md$tumor_stage)
tmp <- mice::mice.impute.polyreg(y = tmp,
                                 ry = !is.na(tmp),
                                 x = model.matrix(~size+
                                                    intclust +grade,
                                                  data = md)[,-1],
                                 wy = array(TRUE, dim = length(tmp)))
md$tumor_stage[is.na(md$tumor_stage)] <- tmp[is.na(md$tumor_stage)]
remove(tmp)
assertthat::assert_that(sum(is.na(md)) == 0)
md$intclust = factor(md$intclust, levels = c("1", "2", "3", "4ER+", "4ER-", "5",
                                             "6", "7", "8", "9", "10"))
md$cohort = factor(md$cohort, levels = c("1", "2", "3", "4", "5"))

#--- Gene matrix preprocess ----- #
gendata <- read_tsv("C:/RFactory/bymetabric_files/metabricdata/brca_metabric/data_expression.txt", col_names = TRUE)
#gendata <- read_tsv("brca_metabric/data_expression.txt", col_names = TRUE)
source("https://bioconductor.org/biocLite.R")
library(Biobase)
###
glimpse(gendata)
gene_names <- gendata %>% dplyr::select(Hugo_Symbol)  #grap gene names
gendata <- gendata %>% dplyr::select(intersect(colnames(gendata), md$patient_id))# get intersection btw clinical and gene values
sample_names <- colnames(gendata) # get sample names
#Center and scale 
gendata <- unname(t(gendata))
colnames(gendata) <- gene_names %>% unlist; rownames(gendata) <- sample_names %>% unlist
preProcgm <-  caret::preProcess(gendata, method = c("center", "scale")) 
metaES <- predict(preProcgm, gendata)
glimpse(metaES)
metaES <- t(metaES)
colnames(metaES) <- sample_names
gendata <- metaES
rownames(gendata) <- gene_names %>% unlist
rm(metaES)
sum(is.na(gendata))
#Convert to expression set
md <- as.data.frame(md %>% filter(patient_id %in% sample_names) %>% slice(match(sample_names, patient_id))) ; row.names(md) <- md$patient_id#requires classical data frame
x <-  as.matrix(gendata) ;colnames(x) <- rownames(md); 
gene_names <- as.data.frame(gene_names);  rownames(gene_names) <- gene_names %>% unlist
brcaES <- Biobase::ExpressionSet(x,
                                 phenoData = as(md, "AnnotatedDataFrame"),
                                 featureData = as(gene_names, "AnnotatedDataFrame"))
assertthat::assert_that(all(md$patient_id == brcaES$patient_id))
rm(list = c("x", "gendata"))
##### Number of genes reduction using Moderated t statistic
library(limma)
fit <- limma::lmFit(brcaES, design = as.integer(as.factor(brcaES$intclust)))
fit <- eBayes(fit)
print(fit)
glimpse(fit$p.value)
test <- exprs(brcaES)[fit$p.value < quantile(fit$p.value, .05),]
gene_names <- as.data.frame(gene_names[fit$p.value < quantile(fit$p.value, .05),])
rownames(gene_names) <- gene_names %>% unlist
brcaES <- Biobase::ExpressionSet(test,
                                 phenoData = as(md, "AnnotatedDataFrame"),
                                 featureData =  as(as.data.frame(gene_names), "AnnotatedDataFrame"))
assertthat::assert_that(all(md$patient_id == brcaES$patient_id))
rm(list = c("test"))
#Imputation using impute Biobase
sum(is.na(exprs(brcaES)))
# require(MSnbase)
# brcaMSN <- MSnbase::as.MSnSet.ExpressionSet(brcaES)
# brcaMSN <- MSnbase::impute(brcaMSN, method = "MinProb")
# Biobase::exprs(brcaES) <- MSnbase::exprs(brcaMSN)
# rm(brcaMSN)
assertthat::assert_that(sum(is.na(Biobase::exprs(brcaES))) == 0)

#save and clean
saveRDS(md,  "Med_Data_Clean.rds")
saveRDS(brcaES, "Gen_Data.rds")
rm(list = ls())

