#source("RandomWalk_ECTD.R")
#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")
#setwd("P:/seojin")

###############
## Functions ##
###############

my.ECTD <- function(i, j, dat)	{
	res <- (dat[i,i]+dat[j,j]-2*dat[i,j])
	return(res)
}

###################
## Normal Sample ##
###################
dat.inv <- read.csv(file.path(sparsematrixdir, "GLASSO_laplacian_invcov_normal.csv"), header = T)
#my.etcd.normal <- outer(1:nrow(dat.inv), 1:ncol(dat.inv), FUN = my.ECTD, vol = dat.normal.v, dat = dat.inv)
my.etcd.normal <- matrix(NA, ncol = ncol(dat.inv), nrow = nrow(dat.inv))
for (i in 1:nrow(dat.inv)) {
	cat(i, "\n")
	for (j in 1:ncol(dat.inv)) {
		my.etcd.normal[i,j] = my.ECTD(i, j, dat = dat.inv)
	}
}
my.etcd.normal <- dat.normal.v * my.etcd.normal
write.csv(my.etcd.normal, file.path(ECTDdir, "ECTD_GLASSO_normal.csv"))

dat.inv <- read.csv(file.path(sparsematrixdir, "MoorePenrose_invcov_normal.csv"), header = T)
#my.etcd.normal <- outer(1:nrow(dat.inv), 1:ncol(dat.inv), FUN = my.ECTD, vol = dat.normal.v, dat = dat.inv)
my.etcd.normal <- matrix(NA, ncol = ncol(dat.inv), nrow = nrow(dat.inv))
for (i in 1:nrow(dat.inv)) {
	cat(i, "\n")
	for (j in 1:ncol(dat.inv)) {
		my.etcd.normal[i,j] = my.ECTD(i, j, dat = dat.inv)
	}
}
my.etcd.normal <- dat.normal.v * my.etcd.normal
write.csv(my.etcd.normal, file.path(ECTDdir, "ECTD_MoorePenrose_normal.csv"))

###################
## Cancer Sample ##
###################
dat.inv <- read.csv(file.path(sparsematrixdir, "GLASSO_laplacian_invcov_cancer.csv"), header = T)
#my.etcd.cancer <- outer(1:nrow(dat.inv), 1:ncol(dat.inv), FUN = my.ECTD, vol = dat.cancer.v, dat = dat.inv)
my.etcd.cancer <- matrix(NA, ncol = ncol(dat.inv), nrow = nrow(dat.inv))
for (i in 1:nrow(dat.inv)) {
	cat(i, "\n")
	for (j in 1:ncol(dat.inv)) {
		my.etcd.cancer[i,j] = my.ECTD(i, j, dat = dat.inv)
	}
}
my.etcd.cancer <- dat.cancer.v * my.etcd.cancer
write.csv(my.etcd.cancer, file.path(ECTDdir, "ECTD_GLASSO_cancer.csv"))

dat.inv <- read.csv(file.path(sparsematrixdir, "MoorePenrose_invcov_cancer.csv"), header = T)
#my.etcd.cancer <- outer(1:nrow(dat.inv), 1:ncol(dat.inv), FUN = my.ECTD, vol = dat.cancer.v, dat = dat.inv)
my.etcd.cancer <- matrix(NA, ncol = ncol(dat.inv), nrow = nrow(dat.inv))
for (i in 1:nrow(dat.inv)) {
	cat(i, "\n")
	for (j in 1:ncol(dat.inv)) {
		my.etcd.cancer[i,j] = my.ECTD(i, j, dat = dat.inv)
	}
}
my.etcd.cancer <- dat.cancer.v * my.etcd.cancer
write.csv(my.etcd.cancer, file.path(ECTDdir, "ECTD_MoorePenrose_cancer.csv"))




