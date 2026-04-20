library(data.table)
library(lme4)
library(glmnet)
library(doRNG)
library(ACAT)
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
dir_output <- "Result_bcac/"

to_plink_id <- function(x) {
  # drop leading "chr"
  x <- sub("^chr", "", x, ignore.case = TRUE)
  # capture chr, pos, ref, alt (allow letters for indels like INS/DEL if present)
  sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+).*", "\\1:\\2:\\3:\\4", x)
}


cell_list <- c('adipose','fibroblast', 'epithelial')

result_total <- vector(mode = "list", length = 3)



# Read gwas_results
gwas_raw <- fread("/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene/bcac_2020_sumstats_for_s-mixcan_hg38.csv")
# gwas_raw[, MarkerName_hg38 := paste0(sub("^chr", "", CHR), ":", POS_hg38, ":", Baseline.Meta, ":", Effect.Meta)]
gwas_raw[, c("V1") := NULL]
gwas_raw <- gwas_raw[
  nchar(Baseline.Meta) == 1 & nchar(Effect.Meta) == 1 &
    toupper(Baseline.Meta) %chin% c("A","C","G","T") &
    toupper(Effect.Meta) %chin% c("A","C","G","T")
]
t=1
for(cell in cell_list){
    # Read MiXcan weights
    # 1) Read input MiXcan weights file
    W <- fread(paste0("subset_weight_mixcan2_", cell, ".csv"), sep = ",", header = TRUE)

    # if predixcan
    # W <- fread(paste0("subset_weight_predixcanlike.csv"), sep = ",", header = TRUE)
    # cell = 'predixcan'


    # 2) Convert var_ID
    W[, MarkerName_hg38 := to_plink_id(xNameMatrix)]
    W[, chr_num := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]
    gene_list = unique(W$ID)
    G = length(gene_list)
    result = data.frame(matrix(ncol = 11, nrow = G))
    colnames(result) <- c('gene_id','chr','type','celltype','MiXcan_snp_num','filtered_snp_num','Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')

    for(g in 1:G){

      result[g,'gene_id'] = gene_list[g]
      result[g,'chr'] = unique(W[which(W$ID == gene_list[g]), 'chr_num'])
      result[g,'type'] = unique(W[which(W$ID == gene_list[g]), 'type'])
      result[g,'celltype'] = cell
      result[g,'MiXcan_snp_num'] = nrow(W[which(W$ID == gene_list[g]), ])

      gene = gene_list[g]
      print(paste("Processing gene:", gene," in ",cell))

      # build filenames from the gene
      bim_path   <- file.path(dir_base, sprintf("%s_snplist_hg38.bim",      gene))
      raw_path   <- file.path(dir_base, sprintf("%s_snplist_hg38_012.raw",  gene))


      # BIM: variant list (6 cols)
      ld_bim <- fread(bim_path, header = FALSE)
      colnames(ld_bim) = c('chr', 'MarkerName_hg38', 'mor', 'POS_hg38', 'A1_bim', 'A2_bim')

      # RAW: genotypes from --recodeA (space-delimited)
      X_ref <- fread(raw_path, sep = " ")
      snp_ids <- ld_bim[[2]]
      X_ref <- X_ref[,7:ncol(X_ref)]
      colnames(X_ref) <- snp_ids

      cat("Loaded:",
          "\n  BIM  :", bim_path, sprintf("(%d variants)", nrow(ld_bim)),
          "\n  RAW  :", raw_path, sprintf("(%d samples x %d cols)", nrow(X_ref), ncol(X_ref)), "\n")

      target_chr <- ld_bim[1,'chr']
      gwas_chr <- gwas_raw[gwas_raw$CHR == paste0("chr", target_chr), ]

      ld_merge = merge(ld_bim, gwas_chr, by='POS_hg38')
      setDT(ld_merge)
      need_flip  <- (ld_merge$A1_bim==ld_merge$Baseline.Meta& ld_merge$A2_bim==ld_merge$Effect.Meta)
      ld_merge[, Beta.meta := ifelse(need_flip, -as.numeric(Beta.meta), as.numeric(Beta.meta))]
      W_s = W[W$ID == gene, ]

      total_input = merge(ld_merge, W_s, by='MarkerName_hg38')
      total_input <- total_input[, c("MarkerName_hg38", "POS_hg38", "chr", "A1_bim", "A2_bim", "Beta.meta", "sdE.meta", "ID","type",
                                     "weight_cell_1", "weight_cell_2",  "stable_id")]
      ref_snp <- total_input$MarkerName_hg38
      X_ref_filtered <- as.matrix(X_ref[, ref_snp, with = FALSE])

      # Build the p x K weight matrix expected by the current SMiXcan package.
      W_mat <- cbind(
        weight_cell_1 = as.numeric(total_input$weight_cell_1),
        weight_cell_2 = as.numeric(total_input$weight_cell_2)
      )

      gwas_results = list()
      gwas_results$Beta = total_input$Beta.meta
      gwas_results$se_Beta = total_input$sdE.meta

      smixcan_result <- SMiXcan_assoc_test_K(
        W = W_mat,
        gwas_results = gwas_results,
        x_g = X_ref_filtered,
        n0 = NULL,
        n1 = NULL,
        family = "gaussian"
      )
      result[g, 'gene_id'] <- gene
      result[g, 'filtered_snp_num'] <- length(ref_snp)
      result[g, c('Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')] <- c(
        smixcan_result$Z_join[1],
        smixcan_result$p_join_vec[1],
        smixcan_result$Z_join[2],
        smixcan_result$p_join_vec[2],
        smixcan_result$p_join
      )
    }
    result_total[[t]] = result
    t=t+1
    write.csv(result,paste0(dir_output, cell,"_bcac_SMiXcan_result.csv"))
    #write.csv(result,paste0(dir_output, 'predixcan',"_bcac_SMiXcan_result.csv"))
  }

merge_result <- rbindlist(result_total)
write.csv(merge_result,paste0(dir_output,"bcac_SMiXcan_result_merged.csv"))

