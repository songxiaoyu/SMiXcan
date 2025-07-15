# --- Stage 1: prepare_mw_gwas_input.R ---

library(dplyr)
library(tidyr)

# Set directories
base_dir <- "/Users/zhusinan/Downloads/S-MiXcan_code_folder"
gwas_dir <- file.path(base_dir, "code_RealData/BCAC/Breast_Cancer_Risk_2020")
output_dir <- "/Users/zhusinan/Downloads/Data2/1000Genome"

# Read MiXcan model weights
mw_input <- read.csv(file.path(base_dir, "code_RealData/MiXcan_recover/MiXcan_model_weights_GTEx_v8_mammary_recover.tsv"),
                     sep = "\t", header = TRUE)

# Parse varID into components
mw_input <- mw_input %>%
  separate(varID, into = c("CHR", "POS", "ref", "alt", "build"), sep = "_")

# Process each chromosome
for (chr in c(12)) {
  cat("Processing chromosome", chr, "\n")
  
  gwas_file <- file.path(gwas_dir, sprintf("chr%d_icogs_hg38.csv", chr))
  if (!file.exists(gwas_file)) {
    cat("  --> File not found, skipping\n")
    next
  }
  
  # Read GWAS summary stats
  gwas_df_raw <- read.csv(gwas_file)
  gwas_df <- gwas_df_raw[, c(1:9, 48)]
  gwas_df$rsid <- sapply(strsplit(gwas_df$SNP.iCOGs, ":"), function(x) x[1])
  gwas_df <- gwas_df[gwas_df$rsid != chr, ]
  colnames(gwas_df)[10] <- "POS"
  
  # Filter model weights for current chromosome
  mw_chr <- mw_input[mw_input$CHR == paste0("chr", chr), ]
  
  # Merge GWAS + model weights
  mw_gwas_input <- merge(gwas_df, mw_chr, by = "POS")
  
  # Save rsid list (optional)
  write.table(data.frame(mw_gwas_input$rsid.x),
              file = file.path(output_dir, sprintf("bcac2020_filtered_chr%d_gwas_id.txt", chr)),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Save full merged object as .rds
  saveRDS(mw_gwas_input, file = file.path(gwas_dir, sprintf("chr%d_mw_gwas_input.rds", chr)))
}
