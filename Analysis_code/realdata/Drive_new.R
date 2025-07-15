
library(data.table)
library(lme4)
library(glmnet)
library(doParallel)
library(doRNG)
library(ACAT)
library(MiXcan)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)

setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder')
lapply(list.files("code_MiXcan_updated", pattern = "\\.R$", full.names = TRUE), source)
lapply(list.files("code_S-MiXcan", pattern = "\\.R$", full.names = TRUE), source)
lapply(list.files("code_Simulation", pattern = "\\.R$", full.names = TRUE)[2], source)


setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData')
# STEP 1 READ DATA


# read Drive pheno data
drive_pheno = read.csv('Drive/oncoarray-drive.pheno.csv')

# read MiXcan weight data
mw_input_raw <- read.csv('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/MiXcan_model_weights_GTEx_v8_mammary_030323.tsv',header=TRUE)
mw_input_raw <- mw_input_raw[, -c(1,2,3,8)]

for(c in 17:21){
#rm(drive_genome)
print("CHR")
print(c)
chr = 22
mw_input <- mw_input_raw[which(mw_input_raw$CHR == sprintf("chr%d", chr)), ]

# read Drive genome data

file_path <- sprintf("Drive/oncoarray_dosages_chr%02d.txt", chr)
coln <- scan(text = readLines(file_path, 1), what = "", quiet = TRUE)[-1]
drive_genome = read.table(file_path)
print("READ GENOME COMPLETED")
colnames(drive_genome)[2:ncol(drive_genome)] <- coln
colnames(drive_genome)[1] <- 'CHR'

#get intersected snp POS
input_pos = mw_input$POS
drive_pos = drive_genome$POS
pos_select = intersect(input_pos, drive_pos) 
print('SNP NUMBER')
print(length(pos_select)) #5260


# total input data
mw_drive_input <- merge(mw_input, drive_genome, by = "POS") 
colnames(mw_drive_input)[16:ncol(mw_drive_input)] = coln[9:length(coln)]

#STEP 2 RUN GWAS
X_genome = t(mw_drive_input[,16:ncol(mw_drive_input)])
D = drive_pheno[drive_pheno$subject_index %in%  rownames(X_genome), "affection_status" ]
# Match rows of drive_pheno to X_genome by subject_index
D <- drive_pheno[match(rownames(X_genome), drive_pheno$subject_index), "affection_status" ]
# run gwas
drive_gwas_result = run_gwas(X_genome, D, 'binomial')
print("GWAS COMPLETED")
drive_gwas_result = data.frame(drive_gwas_result)
# combine gwas results with input
drive_input_total = cbind(drive_gwas_result, mw_drive_input)


file_path <- sprintf("Drive/drive_input_total_chr%02d.rds", chr)

saveRDS(drive_input_total, file_path)
print('SAVE INPUT COMPLETED')

#STEP 3 SPLIT DATA BY GENE
drive_gene <-unique(mw_drive_input$genename) #
drive_split_df <- split(drive_input_total, drive_input_total$genename)


#STEP 4 RUN MIXCAN & S-MIXCAN
print("START RUNNING S-MIXCAN")
drive_result = data.frame(matrix(ncol = 12, nrow = length(drive_gene)))
colnames(drive_result) <- c('p_m_sep_1','p_m_sep_2','p_m_sep','p_m_join_1','p_m_join_2','p_m_join','p_s_sep_1','p_s_sep_2','p_s_sep','p_s_join_1','p_s_join_2','p_s_join')
filtered_list <- list()

for (g in 1:length(drive_gene)) {
  gene = drive_gene[g]
  gene_df <- drive_split_df[[gene]]  # Access each dataframe
  print(paste("Processing gene:", gene))
  # Extract weight columns
  W1 <- gene_df[["weight_cell_1"]]
  W2 <- gene_df[["weight_cell_2"]]
  
  # Filter SNPs with non-zero weights
  selected_snp <- gene_df[which(W1 != 0 | W2 != 0), ]
  # Get genome with snp id
  x_selected = t(selected_snp[,18:ncol(selected_snp)])
  
  # Calculate predicted gene expression value y
  y_predicted <- data.frame(matrix(ncol = 3, nrow = ncol(selected_snp)-17))
  colnames(y_predicted) = c('cell_0', 'cell_1','cell_2')
  # colnames(y_predicted) = c('cell_1','cell_2')
  y_predicted$cell_0 <- rep(1,nrow(x_selected))
  y_predicted$cell_1 <- (as.matrix(x_selected) %*% as.matrix(selected_snp$weight_cell_1))[,1]
  y_predicted$cell_2 <- (as.matrix(x_selected) %*% as.matrix(selected_snp$weight_cell_2))[,1]
  # scale y1, y2
  y_predicted[,2:3] = scale(y_predicted[,2:3])
  y_predicted$subject_index = coln[9:length(coln)]
  
  # run MiXcan association test
  asso_data = merge(y_predicted, drive_pheno, by='subject_index')
  
  drive_mixcan_sep <- MiXcan_association_sep(new_y = asso_data[,c('cell_1','cell_2')], new_cov = data.frame('cov'=rep(0,nrow(x_selected))), new_outcome = as.matrix(asso_data[, 'affection_status']), family  = "binomial")
  drive_mixcan_join <- MiXcan_association_join(new_y = asso_data[,c('cell_1','cell_2')], new_cov = data.frame('cov'=rep(0,nrow(x_selected))), new_outcome = as.matrix(asso_data[, 'affection_status']), family  = "binomial")
  if(is.na(drive_mixcan_join$cell2_p)){
    drive_mixcan_join$cell2_p = drive_mixcan_join$cell1_p
    drive_mixcan_join$p_combined = drive_mixcan_join$cell1_p
  }
  
  drive_result[g, 'p_m_sep_1'] = drive_mixcan_sep$cell1_p
  drive_result[g, 'p_m_sep_2'] = drive_mixcan_sep$cell2_p
  drive_result[g, 'p_m_sep'] = drive_mixcan_sep$p_combined
  drive_result[g, 'p_m_join_1'] = drive_mixcan_join$cell1_p
  drive_result[g, 'p_m_join_2'] = drive_mixcan_join$cell2_p
  drive_result[g, 'p_m_join'] = drive_mixcan_join$p_combined
  
  #run S-MiXcan
  cov_x_g = cov(x_selected)
  
  Y <- as.matrix(y_predicted[,c('cell_0', 'cell_1','cell_2')])
  #Y <- as.matrix(y_predicted[,c('cell_1','cell_2')])
  YtY = t(Y) %*% Y
  gwas_results <- list()
  gwas_results$Beta = selected_snp$Beta
  gwas_results$se_Beta =selected_snp$se_Beta
  family0 = 'binomial'
  
  ft_0 <- glm(affection_status ~0+ cell_0, data=asso_data, family = 'binomial')
  Z_0 <- coef(summary(ft_0))['cell_0', 3]
  
  S_MiXcan_results <-run_S_MiXcan_r(W1, W2, gwas_results, cov_x_g, YtY, family0 = 'binomial', Z_0=Z_0)
  
  drive_result[g, 'p_s_sep_1'] = S_MiXcan_results$p_1_sep
  drive_result[g, 'p_s_sep_2'] = S_MiXcan_results$p_2_sep
  drive_result[g, 'p_s_sep'] = S_MiXcan_results$p_sep
  drive_result[g, 'p_s_join_1'] = S_MiXcan_results$p_1_join
  drive_result[g, 'p_s_join_2'] = S_MiXcan_results$p_2_join
  drive_result[g, 'p_s_join'] = S_MiXcan_results$p_join
}

  rownames(drive_result) <- drive_gene
  file_path <- sprintf("Drive/drive_result_chr%02d.csv", chr)
  write.csv(drive_result , file_path)
  cor(log(drive_result$p_m_join), log(drive_result$p_s_join),use='complete.obs')
}

View(drive_result)
plot(log10(drive_result$p_s_sep), log10(drive_result$p_s_join))
abline(0,1)
plot(drive_result$p_m_sep, drive_result$p_s_sep)
abline(0,1)
plot(drive_result$p_m_sep, drive_result$p_m_join)
abline(0,1)
plot(log10(drive_result$p_m_join), log10(drive_result$p_s_join))
abline(0,1)
cor(log(drive_result$p_m_join), log(drive_result$p_s_join),use='complete.obs')

drive_result_1 <- read.csv("Drive/drive_result_chr05.csv")
plot(drive_result_1$p_m_sep, drive_result_1$p_m_join)
