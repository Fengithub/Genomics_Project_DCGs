#source("Measures_for_Poster.R")

#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")

##########
## Data ##
##########

## Covariance Matrix Summary
cov.normal <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary.csv"), header = T)	
cov.cancer <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary.csv"), header = T)
cov.normal2 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_normal_summary2.csv"), header = T)	
cov.cancer2 <- read.csv(file.path(sparsematrixdir, "GLASSO_step1_cov_cancer_summary2.csv"), header = T)
cov.normal <- as.data.frame(cbind(cov.normal, cov.normal2))
cov.cancer <- as.data.frame(cbind(cov.cancer, cov.cancer2))

## Commuting Time Distance (ECTD) Summary
ctd.normal <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary.csv"), header = T)
ctd.cancer <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary.csv"), header = T)
ctd.normal2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_normal_summary2.csv"), header = T)
ctd.cancer2 <- read.csv(file.path(ECTDdir, "ECTD_MoorePenrose_cancer_summary2.csv"), header = T)
ctd.normal <- as.data.frame(cbind(ctd.normal, ctd.normal2))
ctd.cancer <- as.data.frame(cbind(ctd.cancer, ctd.cancer2))


##############
## Measures ##
##############

gene <- read.csv(file.path(datdir, "tcga_normal_refined_imputed_subdata4000.csv"), header = T)$X

## Covariance Matrix Measures
cov.measure <- data.frame(Gene = gene)
#cov.measure$DensityDiff <- abs(cov.normal$mean - cov.cancer$mean)
#cov.measure$DensityDiffRank <- rank(-cov.measure$DensityDiff)
cov.measure$DensityStat <- abs((cov.normal$mean - cov.cancer$mean)/(cov.normal$mean + cov.cancer$mean))
cov.measure$DensityStatRank <- rank(-cov.measure$DensityStat)
cov.measure$BtwnessStat <- abs(cov.normal$betweenness - cov.cancer$betweenness)
cov.measure$BtwnessStatRank <- rank(-cov.measure$BtwnessStat)
cov.measure$ClonessStat <- abs(cov.normal$closeness - cov.cancer$closeness)
cov.measure$ClonessStatRank <- rank(-cov.measure$ClonessStat)
cov.measure$ClonessStatscale <- abs(cov.normal$closeness-median(cov.normal$closeness) - cov.cancer$closeness + median(cov.cancer$closeness))
cov.measure$ClonessStatscaleRank <- rank(-cov.measure$ClonessStatscale)
#cov.measure$CluCoefDiff <- abs(cov.normal$coef-cov.cancer$coef)
#cov.measure$CluCoefDiffRank <- rank(-cov.measure$CluCoefDiff)
cov.measure$CluCoefStat <- abs(log(cov.normal$coef*cov.normal$mean/(cov.cancer$coef*cov.cancer$mean)))
cov.measure$CluCoefStatRank <- rank(-cov.measure$CluCoefStat)
write.csv(cov.measure, file.path(resultdir, "Cov_GraphicalMeasures.csv"), row.names = F, quote = F)

## Commuting Time Distance (ECTD) Measures
ctd.measure <- data.frame(Gene = gene)
#ctd.measure$DensityDiff <- abs(ctd.normal$mean - ctd.cancer$mean)
#ctd.measure$DensityDiffRank <- rank(-ctd.measure$DensityDiff)
#ctd.measure$DensityRatio <- abs(log(abs(ctd.normal$mean/ctd.cancer$mean)))
#ctd.measure$DensityRatioRank <- rank(-ctd.measure$DensityRatio)
ctd.measure$DensityStat <- abs((ctd.normal$mean - ctd.cancer$mean)/(ctd.normal$mean + ctd.cancer$mean))
ctd.measure$DensityStatRank <- rank(-ctd.measure$DensityStat)
ctd.measure$BtwnessStat <- abs(ctd.normal$betweenness - ctd.cancer$betweenness)
ctd.measure$BtwnessStatRank <- rank(-ctd.measure$BtwnessStat)
ctd.measure$ClonessStat <- abs(ctd.normal$closeness - ctd.cancer$closeness)
ctd.measure$ClonessStatRank <- rank(-ctd.measure$ClonessStat)
ctd.measure$ClonessStatscale <- abs(ctd.normal$closeness-median(ctd.normal$closeness) - ctd.cancer$closeness+median(ctd.cancer$closeness))
ctd.measure$ClonessStatRankscale <- rank(-ctd.measure$ClonessStatscale)
#ctd.measure$CluCoefDiff <- abs(ctd.normal$coef-ctd.cancer$coef)
#ctd.measure$CluCoefDiffRank <- rank(-ctd.measure$CluCoefDiff)
ctd.measure$CluCoefStat <- abs(log(ctd.normal$coef*ctd.normal$mean/(ctd.cancer$coef*ctd.cancer$mean)))
ctd.measure$CluCoefStatRank <- rank(-ctd.measure$CluCoefStat)
write.csv(ctd.measure, file.path(resultdir, "ECTD_GraphicalMeasures.csv"), row.names = F, quote = F)












				