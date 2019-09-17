import numpy as np
import pandas as pd
import pglasso
from sklearn.preprocessing import scale
import csv

AML1 = pd.read_csv("tcga_cancer_refined_imputed_subdata4000.csv", index_col=0,low_memory=False)
AML2 = pd.read_csv("tcga_normal_refined_imputed_subdata4000.csv", index_col=0,low_memory=False)
gene_names = list(AML1.index)
print AML1.shape
print AML2.shape


def load_pathways():
	p0 = open("reactome.txt").readlines()
	p1 = map(lambda x: x.split(":")[1].split(","), p0) # map function means iteration, lambda function do not need return
	p2 = map(lambda x: map(str.strip, x), p1)
	return p2

def map_pathways(gene_names, p2):
	p3 = map(lambda x: map(gene_names.index, x), p2)
	return p3

pathways = load_pathways()
all_genes = []
map(all_genes.extend, pathways)
all_genes = list(set(all_genes))
print len(all_genes)

idxes = np.array([gene_names.index(g) for g in all_genes])
AML1_data = scale(np.array(AML1)[idxes, :].T).T
AML2_data = scale(np.array(AML2)[idxes, :].T).T

pathways = map_pathways(all_genes, pathways)
AML1_S = np.cov(AML1_data)
AML2_S = np.cov(AML2_data)


T1 = pglasso.pglasso(AML1_S, pathways, 0.15)
T2 = pglasso.pglasso(AML2_S, pathways, 0.15)

csvfile = "cancer_lasso.csv"
f = open(csvfile, "wb")
csvwriter = csv.writer(f, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
for i in range(0, T1.shape[0]):
	csvwriter.writerow(T1[i, :])

f.close()

csvfile = "normal_lasso.csv"
f = open(csvfile, "wb")
csvwriter = csv.writer(f, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
for i in range(0, T2.shape[0]):
	csvwriter.writerow(T2[i, :])

f.close()

def get_likelihood(X, T):
	X = np.array(X)
	return np.linalg.slogdet(T)[1] - np.sum(X*T)

print "Train:", get_likelihood(AML1_S, T1)
print "Test:", get_likelihood(AML2_S, T1)
