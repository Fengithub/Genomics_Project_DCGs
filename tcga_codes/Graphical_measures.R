#source("Graphical_measures.R")

#####################
## Set Working Dir ##
#####################
#install.packages("igraph")
library(igraph)
library(PCIT)
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")


##########
## ECTD ##
##########

### Commuting Time
com.normal <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_cancer.csv"), header = T)
#com.normal <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer.csv"), header = T)
#com.normal <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_normal.csv"), header = T)
#com.normal <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal.csv"), header = T)

com.normal <- com.normal[, 2:4001]
com.normal.mean <- apply(com.normal, 2, f <- function(x) mean(x[x>0]))
com.normal.sd <- apply(com.normal, 2, f <- function(x) sd(x[x>0]))
com.normal.min <- apply(com.normal, 2, f <- function(x) min(x[x>0]))
com.normal.max <- apply(com.normal, 2, max)
com.normal.res <- data.frame(mean = com.normal.mean, sd = com.normal.sd, min = com.normal.min, max = com.normal.max)
diag(com.normal) <- 1
com.normal <- 1/(com.normal)
diag(com.normal) <- 0
com.normal.cent.c <- clusteringCoefficient(as.matrix(com.normal))
com.normal.cent.c.global <- clusteringCoefficientPercent(as.matrix(com.normal))
com.normal.cent.c.global
com.normal.res$coef <- com.normal.cent.c

#write.csv(com.normal.res, file.path(ECTDdir, "ECTD_GLASSO_cancer_summary.csv"), row.names = F, quote = F) 
#write.csv(com.normal.res, file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary.csv"), row.names = F, quote = F) 
#write.csv(com.normal.res, file.path(ECTDdir, "ECTD_GLASSO_normal_summary.csv"), row.names = F, quote = F) 
#write.csv(com.normal.res, file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary.csv"), row.names = F, quote = F) 

### graphs
com.normal.res <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_normal_summary.csv"), header = T)
com.normal.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary.csv"), header = T)
com.cancer.res <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_cancer_summary.csv"), header = T)
com.cancer.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary.csv"), header = T)
com.normal2.res <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_normal_summary2.csv"), header = T)
com.normal2.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary2.csv"), header = T)
com.cancer2.res <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_cancer_summary2.csv"), header = T)
com.cancer2.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary2.csv"), header = T)
pairs(com.normal.res)
pairs(com.normal.res2)
pairs(cbind(com.normal.res, com.normal.res2))
pairs(cbind(com.cancer.res, com.cancer.res2))
pairs(cbind(com.normal.res, com.cancer.res))
pairs(cbind(com.normal.res2, com.cancer.res2))
pairs(com.normal2.res)
pairs(com.normal2.res2)
pairs(cbind(com.normal2.res, com.normal2.res2))
pairs(cbind(com.cancer2.res, com.cancer2.res2))
pairs(cbind(com.normal2.res, com.cancer2.res))
pairs(cbind(com.normal2.res2, com.cancer2.res2))

meandiff <- com.normal.res$mean - com.cancer.res$mean
sddiff <- com.normal.res$sd - com.cancer.res$sd
coef.normal <- com.normal.res$coef
coef.cancer <- com.cancer.res$coef
meandiff2 <- com.normal.res2$mean - com.cancer.res2$mean
sddiff2 <- com.normal.res2$sd - com.cancer.res2$sd
coef.normal2 <- com.normal.res2$coef
coef.cancer2 <- com.cancer.res2$coef
pairs(data.frame(meanave = (com.normal.res$mean +  com.cancer.res$mean)/2, meannormal = com.normal.res$mean, meancancer = com.cancer.res$mean, meandiff = meandiff, sddiff = sddiff, coefnormal = coef.normal, coefcancer = coef.cancer,
				meanave2 = (com.normal.res2$mean +  com.cancer.res2$mean)/2, meannormal2 = com.normal.res2$mean, meancancer2 = com.cancer.res2$mean, meandiff2 = meandiff2, sddiff2 = sddiff2, coefnormal2 = coef.normal2, coefcancer2 = coef.cancer2))

par(mfrow = c(2,4))
thres = 2*10^5
hist(com.normal.res$mean[which(com.normal.res$mean < thres)], xlim = c(0, thres), main = "normal: glasso")
hist(com.cancer.res$mean[which(com.cancer.res$mean < thres)], xlim = c(0, thres), main = "cancer: glasso")
hist(com.normal.res2$mean[which(com.normal.res2$mean < thres)], xlim = c(0, thres), main = "normal: moore")
hist(com.cancer.res2$mean[which(com.cancer.res2$mean < thres)], xlim = c(0, thres), main = "cancer: moore")
thres0 = 4*10^5
thres = 9*10^5
hist(com.normal.res$sd[which(com.normal.res$sd < thres)], xlim = c(thres0, thres), main = "normal: glasso")
hist(com.cancer.res$sd[which(com.cancer.res$sd < thres)], xlim = c(thres0, thres), main = "cancer: glasso")
hist(com.normal.res2$sd[which(com.normal.res2$sd < thres)], xlim = c(thres0, thres), main = "normal: moore")
hist(com.cancer.res2$sd[which(com.cancer.res2$sd < thres)], xlim = c(thres0, thres), main = "cancer: moore")

par(mfrow = c(2,4))
thres.clo = 550
thres.clo2 = 3500
hist(com.normal2.res$closeness[com.normal2.res$closeness<thres.clo], main = "normal: glasso")
hist(com.cancer2.res$closeness[com.cancer2.res$closeness<thres.clo], main = "cancer: glasso")
hist(com.normal2.res$betweenness, main = "normal: glasso")
hist(com.cancer2.res$betweenness, main = "cancer: glasso")
hist(com.normal2.res2$closeness[com.normal2.res2$closeness<thres.clo2], main = "normal: Moore")
hist(com.cancer2.res2$closeness[com.cancer2.res2$closeness<thres.clo2], main = "cancer: Moore")
hist(com.normal2.res2$betweenness, main = "normal: Moore")
hist(com.cancer2.res2$betweenness, main = "cancer: Moore")

png("ECTD_all_histogram_plot_minandcoef.png")
par(mfrow = c(2,4))
#thres = 2*10^5
hist(com.normal.res$coef, main = "normal: glasso")
hist(com.cancer.res$coef, main = "cancer: glasso")
hist(com.normal.res2$coef, main = "normal: glasso")
hist(com.cancer.res2$coef, main = "cancer: glasso")
#thres0 = 4*10^5
#thres = 9*10^5
hist(com.normal.res$min, main = "normal: glasso")
hist(com.cancer.res$min, main = "cancer: glasso")
hist(com.normal.res2$min, main = "normal: glasso")
hist(com.cancer.res2$min, main = "cancer: glasso")
dev.off()

res.normal1 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary.csv"), header = T)	 
res.cancer1 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary.csv"), header = T)
com.normal.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary.csv"), header = T)
com.cancer.res2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary.csv"), header = T)
cov.normal <- res.normal1
cov.cancer <- res.cancer1
cmt.normal <- com.normal.res2
cmt.cancer <- com.cancer.res2

png("compareCOVandCMT_basic.png")
res <- data.frame(denstatcov = abs((cov.normal$mean-cov.cancer$mean)/(cov.normal$mean+cov.cancer$mean)),
                  denstatcmt = abs((cmt.normal$mean-cmt.cancer$mean)/(cmt.normal$mean+cmt.cancer$mean)),
                  ccdiffcov = abs(log(cov.normal$coef*cov.normal$mean/(cov.cancer$coef*cov.cancer$mean))),
                  ccdiffcmt = abs(log(cmt.normal$coef*cmt.normal$mean/(cmt.cancer$coef*cmt.cancer$mean))),
                  btwdiff = abs(res.normal$betweenness-res.cancer$betweenness),
                  clodiff = abs(res.normal$closeness - res.cancer$closeness))
dev.off()

################
## Cov Matrix ##
################

cov.dat <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal.csv"), header = T)
cov.dat <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer.csv"), header = T)

diag(cov.dat) <- 0
cov.dat <- abs(cov.dat)
cov.dat.mean <- apply(cov.dat, 2, f <- function(x) mean(x[x>0]))
cov.dat.sd <- apply(cov.dat, 2, f <- function(x) sd(x[x>0]))
cov.dat.min <- apply(cov.dat, 2, f <- function(x) min(x[x>0]))
cov.dat.max <- apply(cov.dat, 2, max)
cov.dat.res <- data.frame(mean = cov.dat.mean, sd = cov.dat.sd, min = cov.dat.min, max = cov.dat.max)
cov.dat.cent.c <- clusteringCoefficient(as.matrix(cov.dat))
cov.dat.cent.c.global <- clusteringCoefficientPercent(as.matrix(cov.dat))
cov.dat.cent.c.global
cov.dat.res$coef <- cov.dat.cent.c

#write.csv(cov.dat.res, file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary.csv"), row.names = F, quote = F) 
write.csv(cov.dat.res, file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary.csv"), row.names = F, quote = F) 

### graphs
com.normal.res <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary.csv"), header = T)
com.cancer.res <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary.csv"), header = T)
pairs(cbind(com.normal.res, com.cancer.res))

############################
## Measurs for Cov Matrix ##
############################
cov.dat <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal.csv"), header = T)
#cov.dat <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer.csv"), header = T)

FNlist <- c("GLASSO_normal", "GLASSO_cancer", "MoorePenrose_normal", "MoorePenrose_cancer")
FNlist <- c("MoorePenrose_normal", "MoorePenrose_cancer")
for (i in 1:2) {
	FN <- file.path(ECTDdir, paste0("ECTD_", FNlist[i], ".csv"))
#cov.dat <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_cancer.csv"), header = T)
#cov.dat <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer.csv"), header = T)
#cov.dat <- read.csv(file.path(ECTDdir, "ECTD_GLASSO_normal.csv"), header = T)
#cov.dat <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal.csv"), header = T)
	cov.dat <- read.csv(FN, header = T)

cov.dat <- cov.dat[,2:4001]
diag(cov.dat) <- 1
cov.dat <- 1/cov.dat
diag(cov.dat) <- 0
cov.dat <- as.matrix(abs(cov.dat))
g <- graph_from_adjacency_matrix(cov.dat, weighted = TRUE)

	my.clo <- closeness(g)
	my.btw <- estimate_betweenness(g, cutoff = -1)
	my.lcc <- transitivity(g, type = "local")
	res.local <- data.frame(closeness = my.clo, betweenness = my.btw, cluscoef = my.lcc)
	write.csv(res.local, file.path(ECTDdir, paste0("ECTD_", FNlist[i], "_summary2.csv")), row.names = F, quote = F)	
	#write.csv(res.local, file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary2.csv"), row.names = F, quote = F)	
	#write.csv(res.local, file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary2.csv"), row.names = F, quote = F)
	
	my.dia <- diameter(g)
	my.rad <- radius(g)
	my.den <- edge_density(g)
	my.gcc <- transitivity(g, type = "global") ## clustering coefficient  0.0001172185 0.0000000000 0.9215903976 1.0000000000
	res.global <- c(my.dia, my.rad, my.den, my.gcc)
	res.global
	write.table(res.global, file.path(sparsematrixdir, paste0("ECTD", FNlist[i],"_globalsummary2.txt")), col.names = F, row.names = F, quote = F)	##
#	write.table(res.global, file.path(sparsematrixdir, "GLASSO_step1_cov_normal_globalsummary2.txt"), col.names = F, row.names = F, quote = F)	## 0.0001172185 0.0000000000 0.9215903976 1.0000000000
#	write.table(res.global, file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_globalsummary2.txt"), col.names = F, row.names = F, quote = F) ## clustering coefficient 0.009033108 0.000000000 0.929287447 1.000000000

}


## Graph
#FNlist <- c("GLASSO_normal", "GLASSO_cancer", "MoorePenrose_normal", "MoorePenrose_cancer")
#read.csv(file.path(sparsematrixdir, paste0(FNlist[i], "_step1_cov_summary2.csv")), header = T)	
res.normal <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary2.csv"), header = T)	
res.cancer <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary2.csv"), header = T)
res.normal1 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary.csv"), header = T)	
res.cancer1 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary.csv"), header = T)
res.normal <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary2.csv"), header = T)	
res.cancer <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary2.csv"), header = T)
res.normal1 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary.csv"), header = T)	
res.cancer1 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary.csv"), header = T)
res <- data.frame(cbind(res.normal, res.cancer))

png("ECTD_cancerVsnormal_summary2_plot.png")
#png("Cov_cancerVsnormal_summary2_plot.png")
pairs(res)
dev.off()

png("ECTD_cancerVsnormal_summary3_plot.png")
#png("Cov_cancerVsnormal_summary3_plot.png")
pairs(data.frame(densitystat = abs(res.normal1$mean-res.cancer1$mean)/(res.normal1$mean+res.cancer1$mean),
      btwstat = abs(res.normal$betweenness-res.cancer$betweenness),
      clostat = abs(res.normal$closeness - res.cancer$closeness),
      coefnormal = res.normal1$coef, coefcancer = res.cancer1$coef,
      coefdiff = abs(res.normal1$coef-res.cancer1$coef)))
dev.off()

png("ECTD_ccstat.png")
#png("Cov_ccstat.png")
plot(1:nrow(res.normal1), abs(log(abs(res.normal1$coef*res.normal1$mean/(res.cancer1$coef*res.cancer1$mean)))))
dev.off()

png("ETCD_cancerVsnormal_summary2_hist.png")
#png("Cov_cancerVsnormal_summary2_hist.png")
par(mfrow = c(2, 3))
hist(res.normal$closeness)
hist(res.normal$betweenness)
hist(res.normal$cluscoef)
hist(res.cancer$closeness)
hist(res.cancer$betweenness)
hist(res.cancer$cluscoef)
dev.off()

clodiff = abs(res.normal$closeness - res.cancer$closeness)
idx3 = which(clodiff>10^-6 & clodiff <= 1.499980*10^{-6});idx3
idx4 = which(clodiff>1.499980*10^-6);idx4
res.normal1 = res.normal1[idx3,]
res.cancer1 = res.cancer1[idx3,]
res.normal = res.normal[idx3,]
res.cancer = res.cancer[idx3,]

png("ETCD_centrality_clogp3_summary2.png", height = 800, width = 800)
#png("Cov_centrality_clogp3_summary2.png", height = 800, width = 800)
#pairs(data.frame(meannormal = res.normal1$mean, meancancer = res.cancer1$mean, meandiff = abs(res.normal1$mean-res.cancer1$mean), sdnormal = res.normal1$sd, sdcancer = res.cancer1$sd, betweenessdiff = abs(res.normal$betweenness - res.cancer$betweenness)))
pairs(data.frame(#sdnormal = res.normal1$sd, sdcancer = res.cancer1$sd, 
                 #densitystat = abs(res.normal1$mean-res.cancer1$mean)/abs(res.normal1$mean+res.cancer1$mean),
                 betweenessdiff = abs(res.normal$betweenness - res.cancer$betweenness),
                 closenessdiff = abs(res.normal$closeness - res.cancer$closeness),
                 closenessnormal = res.normal$closeness, closenesscancer = res.cancer$closeness))
dev.off()

###### Network Measures #####
cov.measure <- read.csv(file.path(resultdir, "Cov_GraphicalMeasures.csv"), header = T)
ctd.measure <- read.csv(file.path(resultdir, "ECTD_GraphicalMeasures.csv"), header = T)

## Comparing All Network Measures
png("All_stats.png", height = 1000, width = 1000)
pairs(data.frame(cbind(cov.measure[,c(2,4,6,8,10)], ctd.measure[,c(2,4,6,8,10)])))
dev.off()

png("All_stats_rank.png", height = 1000, width = 1000)
pairs(data.frame(cbind(cov.measure[,c(3,5,7,9)], ctd.measure[,c(3,5,7,9)])))
dev.off()


par(mfrow = c(2,4))
plot(1:nrow(ctd.cancer), ctd.cancer$closeness)
points(1:nrow(ctd.normal), ctd.normal$closeness, col = "red")





