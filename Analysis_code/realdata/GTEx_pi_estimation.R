library(data.table)
library(xCell)
library(tidyverse)
library(janitor)
library(SMiXcan)
library(readr)
library(dplyr)
library(glmnet)

setwd("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data")

# load GTEx data and clean ----------------------------
gtex_race <- read_csv("gtex_v8_race.csv")
gtex_white <- gtex_race %>% filter(RACE == "White") %>% pull(SUBJID)
#> length(gtex_white)
#[1] 832

cov = data.frame(fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
#> dim(cov)
#[1] 980   4

breast = fread("Breast_Mammary_Tissue.v8.normalized_expression.bed") # Extract the expression data only to save space
#> dim(breast)
#[1] 25849   400

# Add gene name, revise gene id
ensembl38 <-read_csv("ensembl38.txt") %>% clean_names()
ensembl38 = unique(ensembl38[, c("gene_stable_id", "gene_name")])
breast$gene_id2 = matrix(unlist(strsplit(breast$gene_id, '[.]')), ncol = 2, byrow = T)[, 1]
breast2 = merge(ensembl38, breast, by.x = "gene_stable_id", by.y = "gene_id2")
dup = unique(breast2[duplicated(breast2$gene_name), "gene_name"])
breast3 = breast2[-which(breast2$gene_name %in% dup), ]
#> dim(breast3)
# [1] 24919   402

# narrow to 125 white female
exprB = breast3[, colnames(breast3) %in% cov[which(cov$SEX == 2), "SUBJID"]] # We only need female data
exprB <- exprB %>% dplyr::select(intersect(names(exprB), gtex_white)) # only white
rownames(exprB) = breast3$gene_name
dim(exprB)
# [1] 24919   125

# Cell type composition pi estimation -----------------------------------------------------------

# xCell to prior estimate of cell percentage
xCellBScore=xCellAnalysis(expr=exprB)
rownames(xCellBScore)[23] # This is the epithelial cell scores using xCell based on RNAseq data
piPriorB=xCellBScore[23,]
TSNetB_prop <- NULL

#n1 = ncol(exprB)

# get epithelial xCell signature genes
GeneXCell=data.frame(data.table::fread("aran_butte_2017_suppl_3.csv"))
sig=GeneXCell[grep("Epithelial", GeneXCell$Celltype_Source_ID), c(-1, -2)]
sig=unique(as.vector(as.matrix(sig)))
#

GTEx_deidentified_Epithelial_genes <- exprB[rownames(exprB) %in% sig,] %>%
  rownames_to_column(var = "Gene_ID") %>%
  as_tibble()
#set.seed(111)
#sample_name <- sample(1:125)
#names(GTEx_deidentified_Epithelial_genes)[2:126] <- paste0("sample_",sample_name)
#GTEx_deidentified_Epithelial_genes %>% write_tsv("GTEX_Results/GTEx_deidentified_epithelial_genes.tsv")
piPriorB %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(prior = ".") %>%
  mutate(id = paste0("sample_",sample_name)) %>%
  mutate(prior = prior /2.6+0.1) %>%
  write_tsv("GTEX_Results/GTEx_deidentified_prior.tsv")
((piPriorB+ runif(n1, -0.1, 0.1)) /0.8+0.2) %>% summary
#to fix prior datafrmae
prior = as.data.frame(piPriorB/2.6 +0.1)
colnames(prior)[1]=prior

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0793  0.3337  0.5691  0.5298  0.7302  0.9689
#
# pi_estimation using prior
# bootstrap to stablize the algorithm. ### Sinan change to MiXcan pi_estimation
# todo: match id
GTEx_prior = read_tsv('GTEX_Results/GTEx_deidentified_prior.tsv')
GTEx_em = as.matrix(GTEx_deidentified_Epithelial_genes[,-1])
# Make sure column names of expression matrix are sample IDs
head(colnames(GTEx_em))

# Convert 'prior' into a named numeric vector
prior_vec <- setNames(GTEx_prior$prior, GTEx_prior$id)

# Reorder to match the columns of your expression matrix
prior_vec <- prior_vec[colnames(GTEx_em)]

# Double-check alignment
stopifnot(length(prior_vec) == ncol(GTEx_em))
stopifnot(all(names(prior_vec) == colnames(GTEx_em)))


set.seed(123)
# set n_iteration
pis = SMiXcan:pi_estimation(GTEx_em, n_iteration = 100, prior_vec)
#pis = as.data_frame(pis)
pis = as.data.frame(pis)
rownames(pis) = colnames(exprB)
#save pi_GTEx.csv

