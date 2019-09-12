#####################
## Set Working Dir ##
#####################
setwd("~/Box Sync/Genomics_Project/")
source("dirinfo.R")


##########
## Data ##
##########

## File Names
normal.fn <- read.table("tcga_normal_rawdata/FILE_SAMPLE_MAP.txt", header = T)
cancer.fn <- read.table("tcga_cancer_rawdata/FILE_SAMPLE_MAP.txt", header = T)

## Refine Normal Data
fn <- normal.fn$filename[1]
FN <- paste0("tcga_normal_rawdata/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/", fn)
temp.dat <- read.table(FN, header = F, skip = 2)
dat <- data.frame(gene = temp.dat[,1])
for (i in 1:length(normal.fn$filename)) {
	fn <- normal.fn$filename[i]
	FN <- paste0("tcga_normal_rawdata/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/", fn)
	temp.dat <- read.table(FN, header = F, skip = 2)
	if (sum(dim(temp.dat) != c(17814, 2))>0) stop(paste0("error: check data size for file ", fn))
	eval(parse(text = paste0("dat$V", i, "<- temp.dat[,2]"))) 
}
write.table(dat, "tcga_refined_data/tcga_normal_refined.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.csv(dat, "tcga_refined_data/tcga_normal_refined.csv", row.names = F)


## Refine Cancer Data
fn <- cancer.fn$filename[1]
FN <- paste0("tcga_cancer_rawdata/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/", fn)
temp.dat <- read.table(FN, header = F, skip = 2)
dat <- data.frame(gene = temp.dat[,1])
for (i in 1:length(cancer.fn$filename)) {
	fn <- cancer.fn$filename[i]
	FN <- paste0("tcga_cancer_rawdata/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/", fn)
	temp.dat <- read.table(FN, header = F, skip = 2)
	if (sum(dim(temp.dat) != c(17814, 2))>0) stop(paste0("error: check data size for file ", fn))
	eval(parse(text = paste0("dat$V", i, "<- temp.dat[,2]"))) 
}
#dat <- as.matrix(dat)
write.table(dat, "tcga_refined_data/tcga_cancer_refined.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.csv(dat, "tcga_refined_data/tcga_cancer_refined.csv", row.names = F, quote = F)







