#source("GeneNumReduction_tcga.R")

#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")


##########
## Data ##
##########
dat.normal <- read.csv("tcga_refined_data/tcga_normal_refined_imputed.csv")	# 17814    63
dat.cancer <- read.csv("tcga_refined_data/tcga_cancer_refined_imputed.csv")	# 17814   532
my.res <- read.table(file.path(datdir, "tcga_DEGresult.txt"), header = T)


############################
## Reduce the gene number ##
############################

length(which(my.res$BonferPvalue < 0.05)) # 8039
length(which(my.res$BHPvalue < 0.05)) # 12909

## around 8000
sub.genes <- my.res$gene[which(my.res$BonferPvalue < 0.05)]
dat.normal.sub <- dat.normal[match(sub.genes, dat.normal$X),]
dat.cancer.sub <- dat.cancer[match(sub.genes, dat.cancer$X),]
write.csv(dat.normal.sub, file.path(datdir, "tcga_normal_refined_imputed_subdata8039.csv"), row.names = F, quote = F)
write.csv(dat.cancer.sub, file.path(datdir, "tcga_cancer_refined_imputed_subdata8039.csv"), row.names = F, quote = F)

## around 4000
sub.genes <- my.res$gene[which(my.res$BonferPvalue <= sort(my.res$BonferPvalue)[4000])]
dat.normal.sub <- dat.normal[match(sub.genes, dat.normal$X),]
dat.cancer.sub <- dat.cancer[match(sub.genes, dat.cancer$X),]
write.csv(dat.normal.sub, file.path(datdir, "tcga_normal_refined_imputed_subdata4000.csv"), row.names = F, quote = F)
write.csv(dat.cancer.sub, file.path(datdir, "tcga_cancer_refined_imputed_subdata4000.csv"), row.names = F, quote = F)