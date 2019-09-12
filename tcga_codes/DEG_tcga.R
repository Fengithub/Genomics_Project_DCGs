#source("DEG_tcga.R")

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

###############################################
## Find DEG (Differenctially Expressed Gene) ##
###############################################

my.deg.test <- function(i, dat.normal, dat.cancer) {
	gene.normal <- dat.normal[i,-1]
	gene.cancer <- dat.cancer[i,-1]
	gene.test <- t.test(gene.normal, gene.cancer)
	gene.p <- gene.test$p.value
	gene.t <- gene.test$statistic
	gene.res <- c(gene.p, gene.t)
	return(gene.res)
}

if (nrow(dat.normal) == nrow (dat.cancer)) p = nrow(dat.normal)
my.res <- t(sapply(1:p, my.deg.test, dat.normal = dat.normal, dat.cancer = dat.cancer))
my.res <- as.data.frame(my.res)
colnames(my.res) <- c("Pvalue", "Statistic")
my.res$gene <- dat.normal$X
my.res$BHPvalue <- p.adjust(my.res$Pvalue, method = "BH")
my.res$BonferPvalue <- p.adjust(my.res$Pvalue, method = "bonferroni")
write.table(my.res[c("gene", "Statistic", "Pvalue", "BHPvalue", "BonferPvalue")], file.path(datdir, "tcga_DEGresult.txt"), col.names = T, row.names = F, quote = F)