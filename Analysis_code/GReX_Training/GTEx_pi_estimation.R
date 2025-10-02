library(data.table)
library(xCell)
library(tidyverse)
library(janitor)
library(MiXcan)
library(readr)
library(dplyr)
library(glmnet)
library(janitor)
library(tibble)
library(doParallel)
library(dplyr)
# Set working path
setwd("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data")

# Step 1: Load data---------------------
# 1. Load and clean GTEx race data. Select only White people.
gtex_race <- read_csv("gtex_v8_race.csv")
gtex_white <- gtex_race %>% filter(RACE == "White") %>% pull(SUBJID)

# 2. Load data and select breast cancer genes
cov = data.frame(fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
breast = fread("Breast_Mammary_Tissue.v8.normalized_expression.bed") # Extract the expression data only to save space
ensembl38 <-read_csv("ensembl38.txt") %>% clean_names()
ensembl38 = unique(ensembl38[, c("gene_stable_id", "gene_name")])
breast$gene_id2 = matrix(unlist(strsplit(breast$gene_id, '[.]')), ncol = 2, byrow = T)[, 1]
breast2 = merge(ensembl38, breast, by.x = "gene_stable_id", by.y = "gene_id2")
dup = unique(breast2[duplicated(breast2$gene_name), "gene_name"])
breast3 = breast2[-which(breast2$gene_name %in% dup), ]

# 3. Select breast cancer expression level for white female
exprB = breast3[, colnames(breast3) %in% cov[which(cov$SEX == 2), "SUBJID"]] # We only need female data
exprB <- exprB %>% dplyr::select(intersect(names(exprB), gtex_white)) # only white
rownames(exprB) = breast3$gene_name
dim(exprB)
# [1] 24919   125

# Step2: Cell type composition pi estimation -----------------------------------------------------------

# 1. Run xCell to prior estimate of cell percentage
xCellBScore=xCellAnalysis(expr=exprB)
rownames(xCellBScore)[23] # This is the epithelial cell scores using xCell based on RNAseq data
piPriorB=xCellBScore[23,]


# 2. Get epithelial xCell signature genes
GeneXCell=data.frame(data.table::fread("aran_butte_2017_suppl_3.csv"))
sig=GeneXCell[grep("Epithelial", GeneXCell$Celltype_Source_ID), c(-1, -2)]
sig=unique(as.vector(as.matrix(sig)))


GTEx_deidentified_Epithelial_genes <- exprB[rownames(exprB) %in% sig,] %>%
  rownames_to_column(var = "Gene_ID") %>%
  as_tibble()

piPriorB %>%
  as.data.frame() %>%
  as_tibble()

piPriorB <- piPriorB / 2.6 + 0.1

GTEx_em = as.matrix(GTEx_deidentified_Epithelial_genes[,-1])

# Double-check alignment
stopifnot(length(piPriorB) == ncol(GTEx_em))
stopifnot(all(names(piPriorB) == colnames(GTEx_em)))
set.seed(1234)
pis = pi_estimation(GTEx_em, n_iteration = 5, piPriorB)

# output pi score estimate
pis = as.data.frame(pis)
rownames(pis) = colnames(exprB)


