# Run Primo pattern inference on CARDMPRI PWAS S-MiXcan results.

library(data.table)
library(dplyr)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

workspace_dir <- file.path(paper_dir, "Results", "pwas", "cardmpri_workspace")
results_dir <- file.path(workspace_dir, "cardmpri_result")
data_dir <- file.path(paper_dir, "Data")

combined_path <- file.path(results_dir, "cardmpri_result_pwas.csv")
combined <- read.csv(combined_path)

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
write.csv(out, file.path(results_dir, "cardmpri_result_pwas_annotated.csv"), row.names = FALSE)

out_rename <- out %>%
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

write.csv(out_rename, file.path(results_dir, "cardmpri_table_pwas.csv"), row.names = FALSE)
