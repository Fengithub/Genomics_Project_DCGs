#source("ConstructSparseNetwork_from_moore_tcga.R")
#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")


##############
## Packages ##
##############
#install.packages("MASS")
library(MASS)


###################
## Normal Sample ##
###################
dat.net <- read.csv(file.path(sparsematrixdir, "laplacian_normal.csv"), header = T)
aa <- ginv(as.matrix(dat.net))
write.csv(aa, file.path(sparsematrixdir, "MoorePenrose_invcov_normal.csv"), row.names = F)

###################
## Cancer Sample ##
###################

dat.net <- read.csv(file.path(sparsematrixdir, "laplacian_cancer.csv"), header = T)
aa <- ginv(as.matrix(dat.net))
write.csv(aa, file.path(sparsematrixdir, "MoorePenrose_invcov_cancer.csv"), row.names = F)


