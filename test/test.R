library('testthat')
#test-your_function.R
# Load pre-merged data

chr=5

mw_gwas_input <- readRDS(sprintf("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/BCAC/Breast_Cancer_Risk_2020/chr%d_mw_gwas_input.rds", chr))

# Read LD matrix and SNP list
ld_ref <- read.csv(sprintf("/Users/zhusinan/Downloads/Data2/1000Genome/filtered_chr%d_hg38_ld_matrix_r.ld", chr),
                   header = FALSE, sep = "\t")
ld_snp <- read.delim(sprintf("/Users/zhusinan/Downloads/Data2/1000Genome/filtered_chr%d_hg38.bim", chr), header = FALSE)

ref_snp_id <- ld_snp$V2
ld_ref <- as.matrix(ld_ref)
colnames(ld_ref) <- ld_snp$V2
rownames(ld_ref) <- ld_snp$V2
# Read genotype matrix and set proper colnames
X_ref <- read.csv(sprintf("/Users/zhusinan/Downloads/Data2/1000Genome/filtered_chr%d_hg38_012.raw", chr), sep = ' ')
X_ref <- X_ref[, 7:ncol(X_ref)]
colnames(X_ref) <- sapply(strsplit(colnames(X_ref), "_"), `[[`, 1)
X_ref <- as.matrix(X_ref)

# Filter intersecting SNPs
gwas_ref_snps <- intersect(mw_gwas_input$rsid.x, ref_snp_id)
filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$rsid.x %in% gwas_ref_snps, ]
ld_ref <- ld_ref[gwas_ref_snps, gwas_ref_snps]

# Compute allele frequency scaling matrix D and covariance matrix
avg_AF <- aggregate(Freq.Gwas ~ rsid.x, data = filtered_mw_gwas_input, FUN = mean)
d <- (sqrt(2 * avg_AF$Freq.Gwas * (1 - avg_AF$Freq.Gwas)))
cov_ref <- sweep(sweep(ld_ref, 1, d, "*"), 2, d, "*")
cov_ref <- as.matrix(cov_ref)
# cov_ref <- D %*% as.matrix(ld_ref) %*% D
colnames(cov_ref) <- colnames(ld_ref)
rownames(cov_ref) <- colnames(ld_ref)
# Prepare filtered gene list
split_df <- split(filtered_mw_gwas_input, filtered_mw_gwas_input$gene_name)
filtered_list <- list()
for (gene in names(split_df)) {
  gene_df <- split_df[[gene]]
  W1 <- gene_df$weight_cell_1
  W2 <- gene_df$weight_cell_2
  W <- cbind(W1, W2)
  selected_snp <- gene_df[W1 != 0 | W2 != 0, ]
  filtered_list[[gene]] <- list(W = W, selected_snp = selected_snp)
}

# Estimate Z_0
n_case <- 133384
n_control <- 113789 + 18908
D_vec <- c(rep(1, n_case), rep(0, n_control))
y_0 <- rep(1, length(D_vec))
ft_0 <- glm(D_vec ~ 0 + y_0, family = 'binomial')
Z_0 <- coef(summary(ft_0))["y_0", 3]

# Run S-MiXcan
G <- length(filtered_list)
real_result <- data.frame(matrix(ncol = 6, nrow = G))
colnames(real_result) <- c('p_s_sep_1', 'p_s_sep_2', 'p_s_sep', 'p_s_join_1', 'p_s_join_2', 'p_s_join')


test_that("S-MiXcan works correctly", {
  # Define the inputs and expected outputs
  input1 <- 2
  input2 <- 3
  expected_output1 <- 5

  # Test if the function returns the correct sum
  expect_equal(S-MiXcan(input1, input2), expected_output1)

  # Test if the function throws an error for incorrect input
  wrong_input <- "a"
  expect_error(S-MiXcan(wrong_input, input2))
})
