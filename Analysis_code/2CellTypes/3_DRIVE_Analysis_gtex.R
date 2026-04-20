
library(dplyr)
#setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData')
setwd("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data")

# STEP 1 READ DATA-------------
# read MiXcan weight data

mw_input_new <- read.csv("weights_miXcan_full_2025.csv")
colnames(mw_input_new)[2] ='genename'
mw_input_new <- mw_input_new %>%
  mutate(CHR = sub("^(chr[^_]+).*", "\\1", varID))  # chr10, chrX, chrY, chrM

#test<- test[, -c(1,2,3,8)]
#mw_input_old <- mw_input_old[, -c(1,2,3,8)]
# simple
mw_input <- mw_input_new %>%
  filter(weight_cell_1 != 0 | weight_cell_2 != 0)
mw_input <- mw_input_new %>%
  mutate(POS = as.integer(sub("^chr[^_]+_([0-9]+).*", "\\1", varID)))
mw_input <- mw_input[,c("gene_id","genename","CHR","POS","weight_cell_1","weight_cell_2","type")]


# read Drive genome data
coln <- scan(text = readLines('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/Drive/oncoarray_dosages_chr21.txt', 1), what = "", quiet = TRUE)[-1]
drive_genome = read.table('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/DRIVE/oncoarray_dosages_chr21.txt')
colnames(drive_genome)[2:ncol(drive_genome)] <- coln
colnames(drive_genome)[1] <- 'CHR'

# read Drive pheno data
drive_pheno = read.csv('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/DRIVE/oncoarray-drive.pheno.csv')

#get intersected snp POS
input_pos = mw_input$POS
drive_pos = drive_genome$POS
pos_select = intersect(input_pos, drive_pos) #1102/873


# cleaned input data: mw_drive_input
mw_drive_input <- merge(mw_input, drive_genome, by = "POS")
colnames(mw_drive_input)[16:ncol(mw_drive_input)] = coln[9:length(coln)]

#STEP 2 RUN GWAS------------------
X_genome = t(mw_drive_input[,16:ncol(mw_drive_input)])
D = drive_pheno[drive_pheno$subject_index %in%  rownames(X_genome), "affection_status" ]
 # Match rows of drive_pheno to X_genome by subject_index
D <- drive_pheno[match(rownames(X_genome), drive_pheno$subject_index), "affection_status" ]
 # run gwas
drive_gwas_result = run_gwas(X_genome, D, 'binomial')
drive_gwas_result = data.frame(drive_gwas_result)
 # combine gwas results with inpit
drive_input_total = cbind(drive_gwas_result, mw_drive_input)




#STEP 3 RUN MIXCAN & S-MIXCAN--------------
#split DATA BY GENE
drive_gene <-unique(mw_drive_input$genename) #78 90
drive_split_df <- split(drive_input_total, drive_input_total$genename)

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


  #SMiXcan_assoc_test(W1, W2, gwas_results, X_ref_filtered, cov_ref, n0=NULL, n1=NULL, family='gaussian')
  #S_MiXcan_results <-SMiXcan_assoc_test(W1, W2, gwas_results, cov_x_g, YtY, family0 = 'binomial', Z_0=Z_0)

  drive_result[g, 'p_s_sep_1'] = S_MiXcan_results$p_1_sep
  drive_result[g, 'p_s_sep_2'] = S_MiXcan_results$p_2_sep
  drive_result[g, 'p_s_sep'] = S_MiXcan_results$p_sep
  drive_result[g, 'p_s_join_1'] = S_MiXcan_results$p_1_join
  drive_result[g, 'p_s_join_2'] = S_MiXcan_results$p_2_join
  drive_result[g, 'p_s_join'] = S_MiXcan_results$p_join
}

rownames(drive_result) <- drive_gene
write.csv(drive_result, 'drive_result_new.csv')
View(drive_result)
plot(drive_result$p_m_sep, drive_result$p_m_join)
abline(0,1)
plot(drive_result$p_s_sep, drive_result$p_s_join)
abline(0,1)
plot(drive_result$p_m_join, drive_result$p_s_join)
abline(0,1)
