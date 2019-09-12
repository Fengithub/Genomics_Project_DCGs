#source("RandomWalk_tcga.R")

#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")

##########
## Data ##
##########
dat.normal <- read.csv("tcga_sparse_matrix/GLASSO_step1_cov_normal.csv")	# 4000	4000
dat.cancer <- read.csv("tcga_sparse_matrix/GLASSO_step1_cov_cancer.csv")	# 4000	4000

##############
## Analysis ##
##############

## Adjacency Matrix
dat.normal <- abs(dat.normal)
dat.cancer <- abs(dat.cancer)

## Degree Matrix
dat.normal.d <- apply(dat.normal, 2, sum)
dat.cancer.d <- apply(dat.cancer, 2, sum)

## Laplacian Matrix
dat.normal.l <- diag(dat.normal.d) - dat.normal
dat.cancer.l <- diag(dat.cancer.d) - dat.cancer

## Volumn of Graph
dat.normal.v <- sum(dat.normal.d) # 828804.5
dat.cancer.v <- sum(dat.cancer.d) # 732223.2

write.csv(dat.normal.l, "tcga_sparse_matrix/laplacian_normal.csv", row.names = F, quote = F)
write.csv(dat.cancer.l, "tcga_sparse_matrix/laplacian_cancer.csv", row.names = F, quote = F)