library(data.table)

pheno <- fread("/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Data/phenotypes.both_sexes.tsv.bgz")
head(pheno)
dim(pheno)
names(pheno)
# ncase= 360834
# ncontrol = 360