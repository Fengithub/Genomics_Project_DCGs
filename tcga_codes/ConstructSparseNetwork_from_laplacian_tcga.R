#source("ConstructSparseNetwork_from_laplaciantcga.R")
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
dat.net <- read.csv(file.path(sparsematrixdir, "laplacian_normal.csv"), header = T)

# ## Primal Graphicall Lasso
# rho = 0.7 * max(abs(dat.net))
# B <- dpglasso(dat.net, rho = rho, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_normal.csv", row.names = F)
# write.csv(B$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_normal.csv", row.names = F)
# rho.new = rho*0.8
# B.new <- dpglasso(dat.net, X = B$X, invX = B$invX, rho = rho.new, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B.new$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_normal.csv", row.names = F)
# write.csv(B.new$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_normal.csv", row.names = F)
# rm(c(B, B.new))

## Graphical Lasso
a <- glasso(as.matrix(dat.net), rho = 0.1)
write.csv(a$w, file.path(sparsematrixdir, "GLASSO_laplacian_cov_normal.csv"), row.names = F)
write.csv(a$wi, file.path(sparsematrixdir, "GLASSO_laplacian_invcov_normal.csv"), row.names = F)
 
aa <- glasso(as.matrix(dat.net), rho = 0.2, w.init = a$w, wi.init = a$wi)
write.csv(aa$w, file.path(sparsematrixdir, "GLASSO_laplacian_cov_normal.csv"), row.names = F)
write.csv(aa$wi, file.path(sparsematrixdir, "GLASSO_laplacian_invcov_normal.csv"), row.names = F)
rm(c(a, aa))
rm(c("a", "aa"))

###################
## Cancer Sample ##
###################

dat.net <- read.csv(file.path(sparsematrixdir, "laplacian_cancer.csv"), header = T)

# ## Primal Graphicall Lasso
# rho = 0.7 * max(abs(dat.net))
# B <- dpglasso(dat.net, rho = rho, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_cancer.csv", row.names = F)
# write.csv(B$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_cancer.csv", row.names = F)
# rho.new = rho*0.8
# B.new <- dpglasso(dat.net, X = B$X, invX = B$invX, rho = rho.new, outer.Maxiter = 20, outer.tol = 10^6)
# write.csv(B.new$invX, "tcga_sparse_matrix/primalGLASSO_step1_cov_cancer.csv", row.names = F)
# write.csv(B.new$X, "tcga_sparse_matrix/primalGLASSO_step1_invcov_cancer.csv", row.names = F)
# rm(c("B", "B.new"))

## Graphical Lasso
a <- glasso(as.matrix(dat.net), rho = 0.1)

write.csv(a$w, file.path(sparsematrixdir, "GLASSO_laplacian_cov_cancer.csv"), row.names = F)
write.csv(a$wi, file.path(sparsematrixdir, "GLASSO_laplacian_invcov_cancer.csv"), row.names = F)
 
aa <- glasso(as.matrix(dat.net), rho = 0.2, w.init = a$w, wi.init = a$wi)
write.csv(aa$w, file.path(sparsematrixdir, "GLASSO_laplacian_cov_cancer.csv"), row.names = F)
write.csv(aa$wi, file.path(sparsematrixdir, "GLASSO_laplacian_invcov_cancer.csv"), row.names = F)
rm(c(a, aa))
rm(c("a", "aa"))

