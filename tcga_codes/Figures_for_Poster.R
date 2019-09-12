#source("Figures_for_Poster.R")

#####################
## Set Working Dir ##
#####################
#install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
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


## Matrix Measures
cov.measure <- read.csv(file.path(resultdir, "Cov_GraphicalMeasures.csv"), header = T)
ctd.measure <- read.csv(file.path(resultdir, "ECTD_GraphicalMeasures.csv"), header = T)
#Gene DensityStat DensityStatRank BtwnessStat BtwnessStatRank ClonessStat ClonessStatRank ClonessStatscale
#      Gene DensityStat DensityStatRank BtwnessStat BtwnessStatRank ClonessStat ClonessStatRank ClonessStatscale ClonessStatRankscale CluCoefStat CluCoefStatRank
#1  CREB3L1  0.05268447            3477           0            2002    483.3401             982       0.27198879                 2110   0.1630802            3335
#2 C10orf90  0.03660332            3655           0            2002    483.4422             857       0.37416214                 1840   0.1307714            3504

###########
## Graph ##
###########

## Ven Diagram
thres = 20
for (thres in c(10, 20, 30, 40, 50, 100)) {
  venda <- ctd.measure
  venda$DensityStatvenda <- 0
  idx = c(order(ctd.measure[,3])[1:thres])
  venda$DensityStatvenda[idx] <- 1
  venda$BtwnessStatvenda <- 0
  idx = order(cov.measure[,5])[1:thres]
  venda$BtwnessStatvenda[idx] <- 1
  venda$ClonessStatvenda <- 0
  idx = order(ctd.measure[,9])[1:thres]
  venda$ClonessStatvenda[idx] <- 1
  venda$CluCoefStatvenda <- 0
  idx = order(ctd.measure[,11])[1:thres]
  venda$CluCoefStatvenda[idx] <- 1

  myvenda = venda[c("Gene", "DensityStatvenda", "BtwnessStatvenda", "ClonessStatvenda", "CluCoefStatvenda")]
  myvenda$sum = apply(myvenda[,2:5], 1, sum)
  myvenda = myvenda[myvenda$sum>0,]
  myvenda = myvenda[order(-myvenda$sum),]
  apply(myvenda[,2:5], 2, sum)
  dim(myvenda)

  deg.res <- read.table(file.path(datdir, "tcga_DEGresult.txt"), header = T)
  colnames(deg.res) <- c("Gene", "Statistic", "Pvalue", "BHPvalue", "BonferPvalue")
  deg.res$DEGrank <- rank(-abs(deg.res$Statistic))
  myvenda <- merge(myvenda, deg.res, by = "Gene")
  write.csv(myvenda, file.path(resultdir, paste0("VenDiagram_overappedGenes_top",thres,".csv")), row.names = F, quote = F)
}

ven10 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top10.csv"), header = T)
ven20 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top20.csv"), header = T)
ven30 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top30.csv"), header = T)
ven40 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top40.csv"), header = T)
ven50 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top50.csv"), header = T)
ven100 <- read.csv(file.path(resultdir, "VenDiagram_overappedGenes_top100.csv"), header = T)


## Graph of Measures

cov.measure <- read.csv(file.path(resultdir, "Cov_GraphicalMeasures.csv"), header = T)
ctd.measure <- read.csv(file.path(resultdir, "ECTD_GraphicalMeasures.csv"), header = T)

# plot1 Btw vs Den
thre = cov.measure$BtwnessStat[order(-cov.measure$BtwnessStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = cov.measure$BtwnessStat, den = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVSbtwness.png", height = 1000, width = 1000)
p1<-ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Density Stat.", y = "Betweenness Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p1
dev.off()

thre = ctd.measure$DensityStat[order(-ctd.measure$DensityStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = cov.measure$BtwnessStat, btw = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVSbtwness_inv.png", height = 1000, width = 1000)
p1<-ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Density Stat.", x = "Betweenness Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p1
dev.off()

# plot2: Clo vs Den
thre = ctd.measure$ClonessStatscale[order(-ctd.measure$ClonessStatscale)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = ctd.measure$ClonessStatscale, den = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVScloness.png", height = 1000, width = 1000)
p2 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Density Stat.", y = "Closenness Stat. (Scaled)") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p2
dev.off()

thre = ctd.measure$DensityStat[order(-ctd.measure$DensityStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = ctd.measure$ClonessStatscale, btw = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVScloness_inv.png", height = 1000, width = 1000)
p2 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Density Stat.", x = "Closenness Stat. (Scaled)") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p2
dev.off()



# plot3: CC vs Den
thre = ctd.measure$CluCoefStat[order(-ctd.measure$CluCoefStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = ctd.measure$CluCoefStat, den = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVScluscoef.png", height = 1000, width = 1000)
p3 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Density Stat.", y = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p3
dev.off()

thre = ctd.measure$DensityStat[order(-ctd.measure$DensityStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = ctd.measure$CluCoefStat, btw = ctd.measure$DensityStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("densityVScluscoef_inv.png", height = 1000, width = 1000)
p3 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Density Stat.", x = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p3
dev.off()

# plot4: CC vs gene
thre = ctd.measure$CluCoefStat[order(-ctd.measure$CluCoefStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), cc = ctd.measure$CluCoefStat, ccidx = 1:length(ctd.measure$CluCoefStat))
rownames(dat) = dat$gene
dat$Significant <- ifelse(ctd.measure$CluCoefStat >= thre, "Y", "N")
png("ClusteringCoefStatPlot.png", height = 800, width = 1200)
p4<-ggplot(dat, aes(x = ccidx, y = cc)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Clustering Coefficient Stat.", x = "Gene") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 15) +
  geom_text_repel(
    data = subset(dat, cc >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines")
  )
  p4
dev.off() 

# plot5: Btw vs Btw
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = cov.measure$BtwnessStat, den = ctd.measure$DensityStat, btw2 = ctd.measure$BtwnessStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw2 > 10^4, "Y", "N")
png("btwnessVSbtwness.png", height = 1000, width = 1000)
p5<-ggplot(dat, aes(x = btw, y = btw2)) +
  theme(legend.position = "none", axis.text = element_text(size=40)) +
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Betweenness Stat. (Covariance Mat)", y = "Betweenness Stat. (ECTD)") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 30) +
  geom_text_repel(
    data = subset(dat, btw2 > 10^4),
    aes(label = gene),
    size = 25,
    box.padding = unit(0.8, "lines"),
    point.padding = unit(0.3, "lines")
  )
  p5
dev.off() 

# plot6: Density Plot
thre = ctd.measure$DensityStat[order(-ctd.measure$DensityStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), cc = ctd.measure$DensityStat, ccidx = 1:length(ctd.measure$DensityStat))
rownames(dat) = dat$gene
dat$Significant <- ifelse(ctd.measure$DensityStat >= thre, "Y", "N")
png("DensityStatPlot.png", height = 800, width = 1200)
p6<-ggplot(dat, aes(x = ccidx, y = cc)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Density Stat.", x = "Gene") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 15) +
  geom_text_repel(
    data = subset(dat, cc >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines")
  )
  p6
dev.off()

# plot7: Btw Plot
thre = cov.measure$BtwnessStat[order(-cov.measure$BtwnessStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), cc = cov.measure$BtwnessStat, ccidx = 1:length(cov.measure$BtwnessStat))
rownames(dat) = dat$gene
dat$Significant <- ifelse(cov.measure$BtwnessStat >= thre, "Y", "N")
png("BetweennessStatPlot.png", height = 800, width = 1200)
p7<-ggplot(dat, aes(x = ccidx, y = cc)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Betweenness Stat.", x = "Gene") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 15) +
  geom_text_repel(
    data = subset(dat, cc >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines")
  )
  p7
dev.off() 

# plot8: Closeness Plot (Scaled)
thre = ctd.measure$ClonessStatscale[order(-ctd.measure$ClonessStatscale)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), cc = ctd.measure$ClonessStatscale, ccidx = 1:length(ctd.measure$ClonessStatscale))
rownames(dat) = dat$gene
dat$Significant <- ifelse(ctd.measure$ClonessStatscale >= thre, "Y", "N")
png("ClonessStatPlot.png", height = 800, width = 1200)
p8<-ggplot(dat, aes(x = ccidx, y = cc)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Closeness Stat. (Scaled)", x = "Gene") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 15) +
  geom_text_repel(
    data = subset(dat, cc >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines")
  )
p8
dev.off()

# plot9: Closen vs Btwness
thre = ctd.measure$ClonessStatscale[order(-ctd.measure$ClonessStatscale)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = ctd.measure$ClonessStatscale, den = cov.measure$BtwnessStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("btwnessVScloness.png", height = 1000, width = 1000)
p9 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Betweenness Stat.", y = "Closenness Stat. (Scaled)") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p9
dev.off()

thre = cov.measure$BtwnessStat[order(-cov.measure$BtwnessStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = ctd.measure$ClonessStatscale, btw = cov.measure$BtwnessStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("btwnessVScloness_inv.png", height = 1000, width = 1000)
p9 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Betweenness Stat.", x = "Closenness Stat. (Scaled)") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p9
dev.off()

# plot10: CC vs Btwness
thre = ctd.measure$CluCoefStat[order(-ctd.measure$CluCoefStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = ctd.measure$CluCoefStat, den = cov.measure$BtwnessStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("btwnessVScluscoef.png", height = 1000, width = 1000)
p10 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Betweenness Stat.", y = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p10
dev.off()

thre = cov.measure$BtwnessStat[order(-cov.measure$BtwnessStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = ctd.measure$CluCoefStat, btw = cov.measure$BtwnessStat)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("btwnessVScluscoef_inv.png", height = 1000, width = 1000)
p10 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Betweenness Stat.", x = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p10
dev.off()

# plot11: CC vs Closen
thre = ctd.measure$CluCoefStat[order(-ctd.measure$CluCoefStat)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), btw = ctd.measure$CluCoefStat, den = ctd.measure$ClonessStatscale)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("cluscoeffVScloseness.png", height = 1000, width = 1000)
p11 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(x = "Closenness Stat. (Scaled).", y = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p11
dev.off()

thre = ctd.measure$ClonessStatscale[order(-ctd.measure$ClonessStatscale)]
thre = thre[20]
dat <- data.frame(gene = as.character(cov.measure$Gene), den = ctd.measure$CluCoefStat, btw = ctd.measure$ClonessStatscale)
rownames(dat) = dat$gene
dat$Significant <- ifelse(dat$btw >= thre, "Y", "N")
png("cluscoeffVScloseness_inv.png", height = 1000, width = 1000)
p11 <- ggplot(dat, aes(x = den, y = btw)) +
  theme(legend.position = "none") + 
  geom_point(aes(color = Significant), size = 5) +
  labs(y = "Closenness Stat. (Scaled).", x = "Clustering Coefficient Stat.") + 
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 20) +
  geom_text_repel(
    data = subset(dat, btw >= thre),
    aes(label = gene),
    size = 10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
  p11
dev.off()





