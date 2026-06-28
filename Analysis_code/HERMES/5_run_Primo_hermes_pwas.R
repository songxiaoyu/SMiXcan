# Run Primo pattern inference on HERMES PWAS S-MiXcan results.

library(data.table)
library(dplyr)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

workspace_dir <- Sys.getenv(
  "HERMES_WORKSPACE_DIR",
  unset = file.path(paper_dir, "Results", "hermes_pwas", "hermes_workspace")
)
results_dir <- file.path(workspace_dir, "hermes_result")
data_dir <- file.path(paper_dir, "Data")

result_tag <- Sys.getenv("HERMES_RESULT_TAG", unset = "")
result_suffix <- if (nzchar(result_tag)) paste0("_", result_tag) else ""
combined_path <- Sys.getenv(
  "HERMES_RESULT_FILE",
  unset = file.path(results_dir, "hermes_result_pwas.csv")
)
p_join_col <- Sys.getenv("HERMES_P_JOIN_COL", unset = "p_join")
if (!file.exists(combined_path)) {
  stop("Missing HERMES PWAS result file: ", combined_path, call. = FALSE)
}

combined <- read.csv(combined_path)
if (!p_join_col %in% names(combined)) {
  stop("Missing requested HERMES p-value column: ", p_join_col, call. = FALSE)
}
if (p_join_col != "p_join") {
  combined$p_join <- combined[[p_join_col]]
}

ensembl_ref <- read.csv(file.path(data_dir, "ensembl38.txt"))
setDT(combined)
combined[, gene_id_clean := sub("\\..*$", "", gene_id)]
setDT(ensembl_ref)

ensembl_cyto <- unique(
  ensembl_ref[, .(
    ENSG = sub("\\..*$", "", Gene.stable.ID),
    CYTOBAND = Karyotype.band,
    gene_name_ref = Gene.name
  )],
  by = "ENSG"
)

combined <- ensembl_cyto[combined, on = .(ENSG = gene_id_clean)]
combined[, ENSG := NULL]
combined[, gene_name_final := fifelse(
  !is.na(gene_name) & gene_name != "",
  as.character(gene_name),
  as.character(gene_name_ref)
)]
setcolorder(
  combined,
  c(
    "gene_name_final", "gene_id", "CYTOBAND",
    setdiff(names(combined), c("gene_name_final", "gene_id", "CYTOBAND"))
  )
)

combined <- as.data.frame(combined)
primo <- infer_celltype_patterns(
  merged = combined,
  pvals_names = c("p_cardiomyocytes", "p_other"),
  p_join_name = "p_join",
  type_col = "type",
  gene_id_col = "gene_name_final"
)

out <- primo$out
annotated_path <- file.path(results_dir, paste0("hermes_result_pwas", result_suffix, "_annotated.csv"))
write.csv(out, annotated_path, row.names = FALSE)

out_rename <- out %>%
  select(-any_of(c("gene_name", "gene_name_ref", "gene_id_clean"))) %>%
  rename(
    gene_name = gene_name_final,
    model = type,
    Z_cardiomyocyte = Z_cardiomyocytes,
    p_cardiomyocyte = p_cardiomyocytes,
    p_joint = p_join,
    fdr_p_joint = fdr_p_join,
    prob_00 = post_00,
    prob_10 = post_10,
    prob_01 = post_01,
    prob_11 = post_11
  ) %>%
  select(
    gene_name, gene_id, chr, CYTOBAND, model, input_snp_num,
    Z_cardiomyocyte, p_cardiomyocyte, Z_other, p_other,
    p_joint, fdr_p_joint,
    prob_00, prob_10, prob_01, prob_11
  )

table_path <- file.path(results_dir, paste0("hermes_table_pwas", result_suffix, ".csv"))
write.csv(out_rename, table_path, row.names = FALSE)

cat("Input result:", combined_path, "\n")
cat("p_join column:", p_join_col, "\n")
cat("Annotated output:", annotated_path, "\n")
cat("Table output:", table_path, "\n")
