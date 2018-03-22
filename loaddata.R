#download data
library(cgdsr)
library(dplyr)
library(readr)
library(caret)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
list_cancer = getCancerStudies(mycgds)
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[20,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
# Get clinical data for the case list
clinicaldata = getClinicalData(mycgds,mycaselist)
# Get an easier code to read
names(clinicaldata) <- tolower(names(clinicaldata))
clinicaldata <- tibble::rownames_to_column(clinicaldata, var = "patient_id")
clinicaldata$patient_id <- gsub( "\\.","-", clinicaldata$patient_id)
glimpse(clinicaldata)
groupdata <- clinicaldata %>% 
  select(patient_id, intclust,tumor_stage,cohort)
#DFS data
dfmeta <- read_tsv("brca_metabric/Complete_METABRIC_Clinical_Survival_Data__from_METABRIC.txt")
metafeatures <- read.csv("brca_metabric/Complete_METABRIC_Clinical_Features_Data.txt")[,1:13]
glimpse(metafeatures)
dfdata <- cbind(metafeatures, dfmeta)
dfdata <- dfdata %>% select(patient_id,
                            size, grade, lymph_nodes_positive, time, status)
#Join data medical data
md <- left_join(dfdata, groupdata)
glimpse(md)

#save and clean
saveRDS(md , "medicaldfs.rds")
rm(list = ls())

