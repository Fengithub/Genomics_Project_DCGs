#source("impute_tcga_data.R")
#rm(list = ls())

#####################
## Set Working Dir ##
#####################
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")
library(impute)

################
## imputation ##
################
dat.normal <- read.csv(file.path(datdir, "tcga_normal_refined.csv"), header = T, na.strings = "null") 
dat.normal.temp <- dat.normal[,c(2:ncol(dat.normal))]
dat.normal.impute <- impute.knn(as.matrix(dat.normal.temp))$data
dat.normal.impute <- as.data.frame(dat.normal.impute)
rownames(dat.normal.impute) <- dat.normal$gene
write.csv(dat.normal.impute, "tcga_refined_data/tcga_normal_refined_imputed.csv")

dat.cancer <- read.csv(file.path(datdir, "tcga_cancer_refined.csv"), header = T, na.strings = "null") 
dat.cancer.temp <- dat.cancer[,c(2:ncol(dat.cancer))]
dat.cancer.impute <- impute.knn(as.matrix(dat.cancer.temp))$data
dat.cancer.impute <- as.data.frame(dat.cancer.impute)
rownames(dat.cancer.impute) <- dat.cancer$gene
write.csv(dat.cancer.impute, "tcga_refined_data/tcga_cancer_refined_imputed.csv")
