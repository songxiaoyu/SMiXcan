# Run S-MiXcan association for HEARTFAIL using PWAS heart-protein weights.
#
# Inputs from previous steps:
#   2_HEARTFAIL_prepare_data_pwas.R:
#     Results/heartfail_pwas/heartfail_workspace/heartfail_input/*.rds
#   3_1000Genome_keep_eur_plink_heartfail_pwas.sh:
#     Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/*.bim/*.raw
#
# HEARTFAIL is treated as a binary trait by default.
# Set HEARTFAIL_N_CASES and HEARTFAIL_N_CONTROLS before running binomial PWAS.
#
# This active version uses the package-level association function directly.

library(data.table)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

workspace_dir <- Sys.getenv(
  "HEARTFAIL_WORKSPACE_DIR",
  unset = file.path(paper_dir, "Results", "heartfail_pwas", "heartfail_workspace")
)
heartfail_input_dir <- file.path(workspace_dir, "heartfail_input")
ld_input_dir <- file.path(workspace_dir, "heartfail_filtered_id")
result_dir <- file.path(workspace_dir, "heartfail_result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

family <- Sys.getenv("HEARTFAIL_GWAS_FAMILY", unset = "binomial")
if (!family %in% c("binomial", "gaussian")) {
  stop("HEARTFAIL_GWAS_FAMILY must be binomial or gaussian.", call. = FALSE)
}

n1 <- as.integer(Sys.getenv("HEARTFAIL_N_CASES", unset = NA_character_))
n0 <- as.integer(Sys.getenv("HEARTFAIL_N_CONTROLS", unset = NA_character_))
if (family == "binomial" && (is.na(n1) || is.na(n0) || n1 <= 0 || n0 <= 0)) {
  stop(
    "Set HEARTFAIL_N_CASES and HEARTFAIL_N_CONTROLS before running binomial HEARTFAIL PWAS.",
    call. = FALSE
  )
}

parse_bool_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_chr_list <- function(x) {
  as.integer(trimws(unlist(strsplit(x, ","))))
}

chr_list <- parse_chr_list(Sys.getenv("HEARTFAIL_CHR_LIST", unset = paste(1:22, collapse = ",")))
assoc_reg_mode <- Sys.getenv("HEARTFAIL_ASSOC_REG_MODE", unset = "fixed")
if (!assoc_reg_mode %in% c("fixed", "estimate")) {
  stop("HEARTFAIL_ASSOC_REG_MODE must be fixed or estimate.", call. = FALSE)
}

assoc_reg_scale <- as.numeric(Sys.getenv("HEARTFAIL_ASSOC_REG_SCALE", unset = "0.1"))
if (is.na(assoc_reg_scale) || assoc_reg_scale < 0) {
  stop("HEARTFAIL_ASSOC_REG_SCALE must be a non-negative number.", call. = FALSE)
}

gene_log_every <- as.integer(Sys.getenv("HEARTFAIL_GENE_LOG_EVERY", unset = "50"))
max_genes <- as.integer(Sys.getenv("HEARTFAIL_MAX_GENES", unset = "0"))
skip_existing <- parse_bool_env("HEARTFAIL_SKIP_EXISTING", default = FALSE)
write_combined <- parse_bool_env("HEARTFAIL_WRITE_COMBINED", default = TRUE)

cat("Using HEARTFAIL_ASSOC_REG_MODE =", assoc_reg_mode, "\n")
cat("Using HEARTFAIL_ASSOC_REG_SCALE =", assoc_reg_scale, "\n")
cat("Using HEARTFAIL_SKIP_EXISTING =", skip_existing, "\n")
cat("Using HEARTFAIL_WRITE_COMBINED =", write_combined, "\n")

run_gene <- function(gene_id, selected, X_ref, chr) {
  selected_snp_id <- selected$varID
  available <- selected_snp_id[selected_snp_id %in% colnames(X_ref)]
  selected <- selected[selected$varID %in% available, , drop = FALSE]
  if (!nrow(selected)) {
    return(NULL)
  }

  X_ref_filtered <- X_ref[, selected$varID, drop = FALSE]
  W <- cbind(
    selected$weight_cardiomyocytes,
    selected$weight_other
  )
  gwas_results <- list(
    Beta = selected$beta.Gwas,
    se_Beta = selected$SE.Gwas
  )

  fit <- SMiXcan_assoc_test_K(
    W = W,
    gwas_results = gwas_results,
    x_g = X_ref_filtered,
    n0 = n0,
    n1 = n1,
    family = family,
    regularization = assoc_reg_mode,
    reg_scale = assoc_reg_scale
  )

  pre_reg_p_join <- SMiXcan:::safe_ACAT(fit$p_sep)
  scale_col <- if (assoc_reg_mode == "estimate") {
    "p_join_reg_estimate"
  } else {
    paste0("p_join_reg_scale_", gsub("\\.", "p", format(assoc_reg_scale, scientific = FALSE, trim = TRUE)))
  }

  out <- data.frame(
    gene_id = selected$gene_id[1],
    gene_name = selected$gene_name[1],
    chr = chr,
    type = selected$type[1],
    input_snp_num = nrow(selected),
    assoc_reg_mode = assoc_reg_mode,
    assoc_reg_scale_selected = fit$reg_scale_selected,
    assoc_reg_condition = fit$reg_condition,
    Z_cardiomyocytes = fit$Z_join[1],
    p_cardiomyocytes = fit$p_join_vec[1],
    Z_other = fit$Z_join[2],
    p_other = fit$p_join_vec[2],
    p_join = fit$p_join,
    Z_cardiomyocytes_pre_reg = fit$Z_sep[1],
    p_cardiomyocytes_pre_reg = fit$p_sep[1],
    Z_other_pre_reg = fit$Z_sep[2],
    p_other_pre_reg = fit$p_sep[2],
    p_join_pre_reg = pre_reg_p_join,
    Z_cardiomyocytes_reg = fit$Z_join[1],
    p_cardiomyocytes_reg = fit$p_join_vec[1],
    Z_other_reg = fit$Z_join[2],
    p_other_reg = fit$p_join_vec[2],
    p_join_reg = fit$p_join,
    check.names = FALSE
  )
  out[[scale_col]] <- fit$p_join
  out
}

for (chr in chr_list) {
  cat("Processing chromosome", chr, "\n")
  result_path <- file.path(result_dir, sprintf("heartfail_chr%d_result_pwas.csv", chr))
  if (skip_existing && file.exists(result_path)) {
    cat("  --> Existing result found, skipping:", result_path, "\n")
    next
  }

  mw_gwas_input_path <- file.path(
    heartfail_input_dir,
    sprintf("chr%d_mw_gwas_input_heartfail_pwas.rds", chr)
  )
  ld_snp_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_heartfail_pwas.bim", chr))
  x_ref_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_012_heartfail_pwas.raw", chr))

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
    stop("Step 4 input must contain gene_id/gene_name or protein_id/protein_name columns.", call. = FALSE)
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

  filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$varID %in% intersect(ref_snp_id, colnames(X_ref)), ]
  if (nrow(filtered_mw_gwas_input) == 0) {
    cat("  --> No overlap with LD reference SNPs, skipping\n")
    next
  }

  split_df <- split(as.data.frame(filtered_mw_gwas_input), filtered_mw_gwas_input$gene_id)
  if (!is.na(max_genes) && max_genes > 0) {
    split_df <- split_df[seq_len(min(length(split_df), max_genes))]
  }

  result_list <- vector("list", length(split_df))
  for (g in seq_along(split_df)) {
    gene_id <- names(split_df)[g]
    if (gene_log_every > 0 && (g == 1 || g %% gene_log_every == 0 || g == length(split_df))) {
      cat("  Gene", g, "of", length(split_df), ":", gene_id, "\n")
    }
    result_list[[g]] <- run_gene(gene_id, split_df[[g]], X_ref, chr)
  }

  real_result <- rbindlist(result_list, fill = TRUE)
  write.csv(real_result, result_path, row.names = FALSE)
}

if (write_combined) {
  available_results <- file.path(result_dir, sprintf("heartfail_chr%d_result_pwas.csv", chr_list))
  available_results <- available_results[file.exists(available_results)]
  if (length(available_results)) {
    combined <- rbindlist(lapply(available_results, fread), fill = TRUE)
    write.csv(
      combined,
      file.path(result_dir, "heartfail_result_pwas.csv"),
      row.names = FALSE
    )
  }
}
