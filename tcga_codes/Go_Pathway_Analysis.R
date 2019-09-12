# source("GO_Pathway_Analysis.R")

#########################
### Working Directory ###
#########################
setwd("C:\\Users\\Seojin\\Downloads\\")
source("tcga_codes\\dir_info.R")

###############
### Packages ##
###############
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("GOstats", "KEGG.db", "annotate", "genefilter", "org.Hs.eg.db"))
#biocLite(c("GO.db"))
#biocLite(c("mygene"))
library("GO.db")
library("GOstats")
library("annotate")
library("KEGG.db")
library("org.Hs.eg.db")
library("genefilter")
library("mygene")

############
### Data ###
############
biodir = "tcga_Biological_Evaluation"
#dat.all <- read.table("tcga_DEGresult.txt", header = T)
#dat.all <- read.table("tcga_refined_data//tcga_DEGresult.txt", header = T)
#write.table(dat.all$gene, "tcga_genelist1.txt", col.names = F, row.names = F, quote = F)
#dat <- read.delim("tcga_gene_info.txt", header = T)
#dat20 <- read.delim("tcga_top20_summary.txt", header = T)
#dat20 <- dat20[which(dat20$Species=="Homo sapiens"),]
#write.table(dat20, "tcga_top20_info.txt", sep = "\t", col.names = T, row.names = F, quote = F)
#dat <- read.delim(file.path(biodir, "tcga_gene_info.txt"), header = T)
#dat20 <- read.delim(file.path(biodir, "tcga_top20_summary.txt"), header = T)
#dat20 <- dat20[which(dat20$Species=="Homo sapiens"),]
#write.table(dat20, file.path(biodir, "tcga_top20_info.txt"), sep = "\t", col.names = T, row.names = F, quote = F)


## Annotation
geneid <- read.delim("Annotations014850.txt", header = T)
mygeneid.entrez <- geneid$EntrezGeneID[match(dat.all$gene, geneid$GeneSymbol)]
mygeneid <- dat.all$gene
mygeneid <- mygeneid[-match(levels(mygeneid)[1:2], mygeneid)]
mygeneid.entrez2 <- queryMany(mygeneid, scopes="symbol", fields="entrezgene", species="human")

genelist2 <- data.frame(geneSymbol = mygeneid.entrez2$query, geneEntrez = mygeneid.entrez2$entrezgene)
genelist1 <- data.frame(geneSymbol = dat.all$gene, geneEntrez0 = mygeneid.entrez)
genelist <- merge(genelist1, genelist2, by = "geneSymbol")
genelist$Entrez <- mapply(max, genelist$geneEntrez0, genelist$geneEntrez, na.rm = T)
genelist$Entrez[which(genelist$Entrez==-Inf)]=NA
#write.table(genelist[,c(1,4)], file.path(biodir, "tcga_anno_geneid.txt"), col.names = T, row.names = F, quote = F)

############################
## GO & Pathway Analyssis ##
############################
myres <- read.csv(file.path(resultdir, "ECTD_GraphicalMeasures.csv"), header = T)
mygene <- read.table(file.path(biodir, "tcga_anno_geneid.txt"), header = T)
mygene <- mygene[which(is.na(mygene$Entrez)==FALSE),]	# 15386     2
colnames(mygene) <- c("Gene", "Entrez")
myres <- merge(myres, mygene, by = "Gene")
n = 20
geneid = myres$Entrez[c(which(myres$DensityStatRank<=n | myres$BtwnessStatRank<=n | myres$ClonessStatRank<=n | myres$CluCoefStatRank<=n))]

## Pathway Analysis
getGeneSym <- function(pathwayid, hyperGTestResult, annofile) {
	geneid = geneIdsByCategory(hyperGTestResult, catids = pathwayid)[[1]]
	genesym = paste(mget(geneid, annofile), collapse = ", ")
	genesym
}

	geneID_deg = geneid
	geneID_all = unique(mygene$Entrez)
	
	# Get the entrez gene identifiers that are mapped to a GO ID
	x <- org.Hs.egPATH 
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes]) 
	gene_all <- names(xx[which(names(xx) %in% geneID_all)])
	gene_deg <- gene_all[which(gene_all %in% geneID_deg)]

	hgCutoff <- 1
	params <- new("KEGGHyperGParams",
             geneIds= gene_deg,
             universeGeneIds = gene_all,
             annotation = "org.Hs.eg.db",
             pvalueCutoff=hgCutoff,
             testDirection="over")
	KEGG_hgOver <- hyperGTest(params)
	Kegg_hgOver.sum = summary(KEGG_hgOver)
	if (nrow(Kegg_hgOver.sum)>0) {
		Kegg_hgOver.sum$Genes <- sapply(Kegg_hgOver.sum$KEGGID, getGeneSym, hyperGTestResult = KEGG_hgOver, annofile = org.Hs.egSYMBOL)
		resdir = file.path(newdatdir, "sarp_KEGG")
		if (!file.exists(resdir)) dir.create(resdir)
		write.table(Kegg_hgOver.sum, file.path(resultdir, 'tcga_kegg.txt'), row.names = F, col.names = T, sep = "\t", quote = F)
		write.table(Kegg_hgOver.sum$Term, file.path(resultdir, "tcga_kegg_term.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
	}

## GO Analysis

	geneID_deg = geneid
	geneID_all = unique(mygene$Entrez)
	
	# Get the entrez gene identifiers that are mapped to a GO ID
	x <- org.Hs.egGO
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes]) 
	gene_all <- names(xx[which(names(xx) %in% geneID_all)])
	gene_deg <- gene_all[which(gene_all %in% geneID_deg)]

	hgCutoff <- 1
	params <- new("GOHyperGParams",
             geneIds= gene_deg,
             universeGeneIds = gene_all,
             annotation = "org.Hs.eg.db",
			 ontology="BP",
			 pvalueCutoff=hgCutoff,
			 conditional=FALSE,
			 testDirection="over") 		 
	GO_hgOver <- hyperGTest(params)
	GO_hgOver.sum = summary(GO_hgOver)
	paramsCond <- params
	conditional(paramsCond) <- TRUE
	GO_hgCondOver <- hyperGTest(paramsCond)
	GO_hgCondOver.sum=summary(GO_hgCondOver)
	
	if (nrow(GO_hgCondOver.sum)>0) {
		GO_hgCondOver.sum$Genes <- sapply(GO_hgCondOver.sum$GOBPID, getGeneSym, hyperGTestResult = GO_hgCondOver, annofile = org.Hs.egSYMBOL)
		GO_hgCondOver.sum$BHPvalue <- p.adjust(GO_hgCondOver.sum$Pvalue, method = "BH")
		GO_hgCondOver.sum <- GO_hgCondOver.sum[,c("GOBPID", "Pvalue", "BHPvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term", "Genes")]
		GO_hgCondOver.sum2 <- GO_hgCondOver.sum[which(GO_hgCondOver.sum$Pvalue<0.05 & GO_hgCondOver.sum$Count > 1), ] ##remove raw p value exceed 0.05, or the count is only 1
		resdir = file.path(newdatdir, "sarp_GO")
		if (!file.exists(resdir)) dir.create(resdir)
		write.table(GO_hgCondOver.sum, file.path(resultdir, "tcga_go.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
		write.table(GO_hgCondOver.sum$Term, file.path(resultdir, "tcga_go_term.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
		write.table(GO_hgCondOver.sum2, file.path(resultdir, "tcga_go2.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
		}


