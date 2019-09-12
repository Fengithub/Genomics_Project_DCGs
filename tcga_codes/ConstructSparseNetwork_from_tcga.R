#source("ConstructSparseNetwork_from_tcga.R")
#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")


##############
## Packages ##
##############
#install.packages("glasso")
library(dpglasso)
library(glasso)


###################
## Normal Sample ##
###################
dat.normal <- read.csv(file.path(datdir, "tcga_normal_refined_imputed_subdata4000.csv"), header = T)
dat.normal.geneid <- dat.normal$gene
dat.normal <- t(dat.normal[,c(2:ncol(dat.normal))])

# ## Primal Graphicall Lasso
# dat.net <- cov(dat.normal, use = "pairwise.complete.obs")
# #dat.net <- dat.net[row(dat.net)>col(dat.net)]
# rm("dat.normal")
# rho = 0.7 * max(abs(dat.net))
# B <- dpglasso(dat.net, rho = rho, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_normal.csv", row.names = F)
# write.csv(B$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_normal.csv", row.names = F)
 
rho.new = rho*0.8
B.new <- dpglasso(dat.net, X = B$X, invX = B$invX, rho = rho.new, outer.Maxiter = 20, outer.tol = 10^6)
write.csv(B.new$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_normal.csv", row.names = F)
write.csv(B.new$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_normal.csv", row.names = F)
rm(c("B", "B.new"))

## Graphical Lasso
a <- glasso(dat.net, rho = 0.1)
write.csv(a$w, "tcga_sparse_matrix/GLASSO_step1_cov_normal.csv", row.names = F)
write.csv(a$wi, "tcga_sparse_matrix/GLASSO_step1_invcov_normal.csv", row.names = F)
 
aa <- glasso(dat.net, rho = 0.2, w.init = a$w, wi.init = a$wi)
write.csv(aa$w, "tcga_sparse_matrix/GLASSO_step1_cov_normal.csv", row.names = F)
write.csv(aa$wi, "tcga_sparse_matrix/GLASSO_step1_invcov_normal.csv", row.names = F)
rm(c("a", "aa"))


###################
## Cancer Sample ##
###################
dat.cancer <- read.csv(file.path(datdir, "tcga_cancer_refined_imputed_subdata4000.csv"), header = T)
dat.cancer.geneid <- dat.cancer$gene
dat.cancer <- t(dat.cancer[,c(2:ncol(dat.cancer))])

# ## Primal Graphicall Lasso
# dat.net <- cov(dat.cancer, use = "pairwise.complete.obs")
# rm("dat.cancer")
## dat.net <- dat.net[row(dat.net)>col(dat.net)]
# rho = 0.7 * max(abs(dat.net))
# B <- dpglasso(dat.net, rho = rho, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_cancer.csv", row.names = F)
# write.csv(B$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_cancer.csv", row.names = F)

rho.new = rho*0.8
B.new <- dpglasso(dat.net, X = B$X, invX = B$invX, rho = rho.new, outer.Maxiter = 20, outer.tol = 10^6)
write.csv(B.new$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_cancer.csv", row.names = F)
write.csv(B.new$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_cancer.csv", row.names = F)
rm(c("B", "B.new"))

## Graphical Lasso
a <- glasso(dat.net, rho = 0.1)
write.csv(a$w, "tcga_sparse_matrix/GLASSO_step1_cov_cancer.csv", row.names = F)
write.csv(a$wi, "tcga_sparse_matrix/GLASSO_step1_invcov_cancer.csv", row.names = F)
 
aa <- glasso(dat.net, rho = 0.2, w.init = a$w, wi.init = a$wi)
write.csv(aa$w, "tcga_sparse_matrix/GLASSO_step1_cov_cancer.csv", row.names = F)
write.csv(aa$wi, "tcga_sparse_matrix/GLASSO_step1_invcov_cancer.csv", row.names = F)
rm(c("a", "aa"))


