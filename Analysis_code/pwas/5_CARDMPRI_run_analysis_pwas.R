# Run S-MiXcan association for CARDMPRI using PWAS weights.
#
# Inputs from previous steps:
#   3_CARDMPRI_prepare_data_pwas.R:
#     Results/pwas/cardmpri_workspace/cardmpri_input/*.rds
#   4_1000Genome_keep_eur_plink_pwas.sh:
#     Results/pwas/cardmpri_workspace/cardmpri_filtered_id/*.bim/*.raw
#
# CARDMPRI is treated as a binary trait by default.
# Default sample counts from Data/phenotypes.both_sexes.tsv.bgz:
#   cases    = 360834
#   controls = 360
# Override with PWAS_N_CASES and PWAS_N_CONTROLS if these counts change.
# If you want gaussian smoke-test mode, set PWAS_GWAS_FAMILY=gaussian.

library(data.table)
library(dplyr)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

workspace_dir <- file.path(paper_dir, "Results", "pwas", "cardmpri_workspace")
cardmpri_input_dir <- file.path(workspace_dir, "cardmpri_input")
ld_input_dir <- file.path(workspace_dir, "cardmpri_filtered_id")
result_dir <- file.path(workspace_dir, "cardmpri_result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

family <- Sys.getenv("PWAS_GWAS_FAMILY", unset = "binomial")
if (!family %in% c("binomial", "gaussian")) {
  stop("PWAS_GWAS_FAMILY must be binomial or gaussian.", call. = FALSE)
}

n1 <- as.integer(Sys.getenv("PWAS_N_CASES", unset = "360834"))
n0 <- as.integer(Sys.getenv("PWAS_N_CONTROLS", unset = "360"))
if (family == "binomial" && (is.na(n1) || is.na(n0) || n1 <= 0 || n0 <= 0)) {
  stop(
    "Set PWAS_N_CASES and PWAS_N_CONTROLS before running binomial CARDMPRI PWAS.",
    call. = FALSE
  )
}

chr_list <- as.integer(strsplit(Sys.getenv("PWAS_CHR_LIST", unset = paste(1:22, collapse = ",")), ",")[[1]])

for (chr in chr_list) {
  cat("Processing chromosome", chr, "\n")

  mw_gwas_input_path <- file.path(
    cardmpri_input_dir,
    sprintf("chr%d_mw_gwas_input_cardmpri_pwas.rds", chr)
  )
  ld_snp_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_pwas.bim", chr))
  x_ref_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_012_pwas.raw", chr))

  if (!file.exists(mw_gwas_input_path) || !file.exists(ld_snp_path) || !file.exists(x_ref_path)) {
    cat("  --> Missing input files, skipping\n")
    next
  }

  mw_gwas_input <- readRDS(mw_gwas_input_path)
  if (!"gene_id" %in% names(mw_gwas_input) && "protein_id" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_id <- mw_gwas_input$protein_id
  }
  if (!"gene_name" %in% names(mw_gwas_input) && "protein_name" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_name <- mw_gwas_input$protein_name
  }
  if (!all(c("gene_id", "gene_name") %in% names(mw_gwas_input))) {
    stop("Step 5 input must contain gene_id/gene_name or protein_id/protein_name columns.", call. = FALSE)
  }
  if (nrow(mw_gwas_input) == 0) {
    cat("  --> Empty merged weight/GWAS input, skipping\n")
    next
  }

  ld_snp <- fread(ld_snp_path, header = FALSE)
  ref_snp_id <- ld_snp$V2

  x_ref_dt <- fread(x_ref_path)
  if (ncol(x_ref_dt) <= 6) {
    cat("  --> Empty LD reference raw file, skipping\n")
    next
  }
  X_ref <- as.matrix(x_ref_dt[, 7:ncol(x_ref_dt)])
  colnames(X_ref) <- sub("_[^_]+$", "", colnames(X_ref))

  gwas_ref_snps <- intersect(mw_gwas_input$varID, ref_snp_id)
  filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$varID %in% gwas_ref_snps, ]
  if (nrow(filtered_mw_gwas_input) == 0) {
    cat("  --> No overlap with LD reference SNPs, skipping\n")
    next
  }

  split_df <- split(as.data.frame(filtered_mw_gwas_input), filtered_mw_gwas_input$gene_id)
  filtered_list <- list()
  for (gene_id in names(split_df)) {
    protein_df <- split_df[[gene_id]]
    W <- cbind(
      protein_df$weight_cardiomyocytes,
      protein_df$weight_other
    )
    filtered_list[[gene_id]] <- list(W = W, selected_snp = protein_df)
  }

  real_result <- data.frame(matrix(ncol = 10, nrow = length(filtered_list)))
  colnames(real_result) <- c(
    "gene_id", "gene_name", "chr", "type", "input_snp_num",
    "Z_cardiomyocytes", "p_cardiomyocytes", "Z_other", "p_other", "p_join"
  )
  real_result$chr <- chr

  for (g in seq_along(filtered_list)) {
    gene_id <- names(filtered_list)[g]
    cat("  Gene:", gene_id, "\n")

    selected <- filtered_list[[gene_id]]$selected_snp
    selected_snp_id <- selected$varID
    X_ref_filtered <- X_ref[, selected_snp_id, drop = FALSE]

    gwas_results <- list(
      Beta = selected$beta.Gwas,
      se_Beta = selected$SE.Gwas
    )

    fit <- SMiXcan_assoc_test_K(
      W = filtered_list[[gene_id]]$W,
      gwas_results = gwas_results,
      x_g = X_ref_filtered,
      n0 = n0,
      n1 = n1,
      family = family
    )

    real_result[g, c("gene_id", "gene_name", "chr", "type", "input_snp_num")] <- c(
      selected$gene_id[1],
      selected$gene_name[1],
      selected$CHR[1],
      selected$type[1],
      nrow(filtered_list[[gene_id]]$W)
    )
    real_result[g, c("Z_cardiomyocytes", "p_cardiomyocytes", "Z_other", "p_other", "p_join")] <- c(
      fit$Z_join[1],
      fit$p_join_vec[1],
      fit$Z_join[2],
      fit$p_join_vec[2],
      fit$p_join
    )
  }

  result_path <- file.path(result_dir, sprintf("cardmpri_chr%d_result_pwas.csv", chr))
  write.csv(real_result, result_path, row.names = FALSE)
}

available_results <- file.path(result_dir, sprintf("cardmpri_chr%d_result_pwas.csv", chr_list))
available_results <- available_results[file.exists(available_results)]
if (length(available_results)) {
  combined <- do.call(rbind, lapply(available_results, read.csv))
  write.csv(
    combined,
    file.path(result_dir, "cardmpri_result_pwas.csv"),
    row.names = FALSE
  )
}
