library(data.table)
library(lme4)
library(glmnet)
library(doRNG)
library(ACAT)
library(MiXcan)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)
library(SMiXcan)


setwd('/Users/zhusinan/Downloads/adriana/')
dir_base <- "/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene"
# Filtered reference genome and ld matrix path
dir_base <- 'Filtered_Ref/'
# Output path
dir_output <- "Result_PrediXcan/"

to_plink_id <- function(x) {
  # drop leading "chr"
  x <- sub("^chr", "", x, ignore.case = TRUE)
  # capture chr, pos, ref, alt (allow letters for indels like INS/DEL if present)
  sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+).*", "\\1:\\2:\\3:\\4", x)
}

# ,
pheno_list <- c('DA', 'NDA', 'PMD')
#cell_list <- c('adipose','fibroblast', 'epithelial')

result_total <- vector(mode = "list", length = 3)
t=1

# Read MiXcan weights
# 1) Read input MiXcan weights file
W <- fread(paste0("subset_weight_predixcanlike.csv"), sep = ",", header = TRUE)
W[, MarkerName_hg38 := to_plink_id(xNameMatrix)]
W[, chr_num := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]
gene_list = unique(W$ID)
G = length(gene_list)


for(pheno in pheno_list){
  # Read gwas_results
  gwas_raw <- fread(paste0("MD_results_2021/", pheno, "_MetaResultswInfo_20210609_hg38.txt"))
  gwas_filtered <- gwas_raw[
    nchar(Allele1) == 1 & nchar(Allele2) == 1 &
      toupper(Allele1) %chin% c("A","C","G","T") &
      toupper(Allele2) %chin% c("A","C","G","T")
  ]

  result = data.frame(matrix(ncol = 11, nrow = G))
  colnames(result) <- c('gene_id','chr','type','phenotype','PrediXcan_snp_num','filtered_snp_num','Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')

  for(g in 1:G){

    result[g,'gene_id'] = gene_list[g]
    result[g,'chr'] = unique(W[which(W$ID == gene_list[g]), 'chr_num'])
    result[g,'type'] = unique(W[which(W$ID == gene_list[g]), 'type'])
    result[g,'phenotype'] = pheno
    result[g,'PrediXcan_snp_num'] = nrow(W[which(W$ID == gene_list[g]), ])

    gene = gene_list[g]
    print(paste("Processing gene:", gene, " for ", pheno))

    # build filenames from the gene
    bim_path   <- file.path(dir_base, sprintf("%s_snplist_hg38.bim",      gene))
    raw_path   <- file.path(dir_base, sprintf("%s_snplist_hg38_012.raw",  gene))

    # --- read files ---

    # BIM: variant list (6 cols)
    ld_bim <- fread(bim_path, header = FALSE)
    colnames(ld_bim) = c('CHR', 'MarkerName_hg38', 'mor', 'POS_hg38', 'A1_bim', 'A2_bim')

    # RAW: genotypes from --recodeA (space-delimited)
    X_ref <- fread(raw_path, sep = " ")

    # (optional) set LD dimnames from BIM col2 so you can align by ID
    snp_ids <- ld_bim[[2]]
    stopifnot(nrow(ld_bim) == length(snp_ids))

    X_ref <- X_ref[,7:ncol(X_ref)]
    colnames(X_ref) <- snp_ids

    cat("Loaded:",
        "\n  BIM  :", bim_path, sprintf("(%d variants)", nrow(ld_bim)),
        "\n  RAW  :", raw_path, sprintf("(%d samples x %d cols)", nrow(X_ref), ncol(X_ref)), "\n")

    target_chr <- ld_bim[1,'CHR']
    gwas_chr <- gwas_filtered[gwas_filtered$chr == paste0("chr", target_chr), ]

    ld_merge = merge(ld_bim, gwas_chr, by='POS_hg38')
    ld_merge$Allele1 <- toupper(ld_merge$Allele1)
    ld_merge$Allele2 <- toupper(ld_merge$Allele2)
    setDT(ld_merge)
    need_flip  <- (ld_merge$A1_bim==ld_merge$Allele2 & ld_merge$A2_bim==ld_merge$Allele1)
    ld_merge[, Zscore := ifelse(need_flip, -as.numeric(Zscore), as.numeric(Zscore))]

    W_s = W[W$ID == gene, ]

    total_input = merge(ld_merge, W_s, by='MarkerName_hg38')
    total_input <- total_input[, c("MarkerName_hg38", "POS_hg38", "chr", "A1_bim", "A2_bim", "Zscore", "P-value", "ID","type",
                               "weight_cell_1", "weight_cell_2", "freq", "stable_id")]
    ref_snp <- total_input$MarkerName_hg38
    X_ref_filtered <- as.matrix(X_ref[, ref_snp, with = FALSE])

    # Extract weight columns
    W1 <- as.matrix(total_input$weight_cell_1)
    W2 <- as.matrix(total_input$weight_cell_2)

    gwas_results = list()
    gwas_results$Beta = total_input$Zscore
    gwas_results$se_Beta = rep(1,length(total_input$Zscore))

    SMiXcan_result <- SMiXcan_assoc_test(W1, W2, gwas_results, X_ref_filtered, n0=NULL, n1=NULL, family='gaussian')
    result[g, 'gene_id'] <- gene
    result[g, 'filtered_snp_num'] <- length(ref_snp)
    result[g, c('Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')] <- SMiXcan_result
  }
    result_total[[t]] = result
    t=t+1
  write.csv(result,paste0(dir_output, pheno,"_SMiXcan_result_predixcan_weights.csv"))
}

merge_result <- rbindlist(result_total)
write.csv(merge_result,paste0(dir_output,"SMiXcan_result_merged_predixcan_weights.csv"))
