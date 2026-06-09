# PWAS CARDMPRI Pipeline

This folder contains the PWAS workflow for GTEx heart proteomics using two cell
components:

1. `Cardiomyocytes`
2. `Other`, defined as Fibroblasts + Endothelial + Smooth muscle/pericyte + Immune

The CARDMPRI GWAS input is the hg38-lifted UK Biobank file:

```text
Heart/Data/I9_CARDMPRI.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
```

The old MI GWAS liftover script is kept separately as `liftover_mi_gwas.py`; it
is not part of the CARDMPRI PWAS workflow.

## Scripts

```text
1_liftover_cardmpri_gwas.py
2_train_model_pwas.R
3_CARDMPRI_prepare_data_pwas.R
4_1000Genome_keep_eur_plink_pwas.sh
5_CARDMPRI_run_analysis_pwas.R
6_run_Primo_pwas.R
liftover_mi_gwas.py
```

## Step 1. Liftover CARDMPRI GWAS

This step lifts the original CARDMPRI GWAS to hg38 and adds rsID annotation.
The current downstream pipeline uses the already-created output:

```text
Heart/Data/I9_CARDMPRI.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
Heart/Data/I9_CARDMPRI.gwas.imputed_v3.both_sexes_hg38_rsid_chr_pos.txt
```

Run only if the hg38 CARDMPRI GWAS needs to be regenerated.

## Step 2. Train PWAS Weights

Script:

```bash
Rscript Analysis_code/pwas/2_train_model_pwas.R
```

Main inputs:

```text
Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData
Heart/GTEx_Pi_Estimate/BayesDeBulk_pi.tsv
New generated files/covariate_EA_with_age.txt
New generated files/codes/pruned_by_chr/chr*_dosage_nomiss_LDpruned_500kb_1_r2_0.8.raw
New generated files/codes/pruned_by_chr/pruned_variants_500kb_1_r2_0.8.annotation.tsv
```

The protein RData contains object `Imputed_protein`, a data frame with:

```text
gene_name, gene_id, chr, start, end, GTEx sample columns...
```

The output uses `gene_id` and `gene_name`, not `protein_id`, because the protein
abundance features are indexed by Ensembl gene IDs.

Outputs:

```text
Results/pwas/training_model_weights/weights_pwas_cardiomyocytes_other.csv
Results/pwas/training_model_weights/weights_pwas_cardiomyocytes_other_skipped.csv
```

Per-chromosome outputs are also supported:

```bash
PWAS_CHR_FILTER=1 PWAS_OUTPUT_SUFFIX=_chr1 Rscript Analysis_code/pwas/2_train_model_pwas.R
```

Useful variables:

```text
PWAS_CHR_FILTER       Restrict training to one chromosome.
PWAS_OUTPUT_SUFFIX    Add a suffix to output file names.
PWAS_MAX_PROTEINS     Limit number of genes for a test run.
PWAS_GENO_RAW_DIR     Override pruned dosage raw directory.
PWAS_RSID_ANNOT_FILE  Override pruned rsID/position/allele annotation file.
```

Combined weights format:

```text
gene_id
gene_name
varID
chr
pos
ref_allele
eff_allele
dosed_allele
weight_cardiomyocytes
weight_other
type
```

## Step 3. Prepare CARDMPRI GWAS + Weight Inputs

Script:

```bash
Rscript Analysis_code/pwas/3_CARDMPRI_prepare_data_pwas.R
```

This step merges trained PWAS weights with CARDMPRI GWAS summary statistics by
variant ID. It handles both forward and reversed allele orientation. If the GWAS
alleles are reversed relative to the model weight `varID`, the GWAS beta is
multiplied by `-1`.

Inputs:

```text
Results/pwas/training_model_weights/weights_pwas_cardiomyocytes_other.csv
Heart/Data/I9_CARDMPRI.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
```

Outputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_input/chr*_mw_gwas_input_cardmpri_pwas.rds
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/cardmpri_filtered_chr*_gwas_id_pwas.txt
```

Run a single chromosome:

```bash
PWAS_CHR_LIST=1 Rscript Analysis_code/pwas/3_CARDMPRI_prepare_data_pwas.R
```

Use a specific weights file:

```bash
PWAS_CHR_LIST=1 \
PWAS_WEIGHTS_FILE="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/pwas/training_model_weights/weights_pwas_cardiomyocytes_other_chr1.csv" \
Rscript Analysis_code/pwas/3_CARDMPRI_prepare_data_pwas.R
```

## Step 4. Build 1000Genome EUR LD Reference

Script:

```bash
bash Analysis_code/pwas/4_1000Genome_keep_eur_plink_pwas.sh
```

This step uses the SNP lists from step 3 to extract matching SNPs from 1000Genome
hg38 pfiles, keeps EUR samples, creates PLINK bed/bim/fam files, and exports an
additive `.raw` dosage matrix for LD estimation.

Inputs:

```text
Data/plink_snplist_by_gene/chr*_hg38.pgen
Data/plink_snplist_by_gene/chr*_hg38.pvar
Data/plink_snplist_by_gene/chr*_hg38.psam
Data/1000Genome/eur_ids.txt
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/cardmpri_filtered_chr*_gwas_id_pwas.txt
```

Outputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_pwas.bed
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_pwas.bim
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_pwas.fam
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_012_pwas.raw
```

PLINK tools are installed in conda environment `plink2`:

```bash
PWAS_PLINK2_ENV=plink2 PWAS_PLINK_ENV=plink2 PWAS_CHR_LIST=1 \
bash Analysis_code/pwas/4_1000Genome_keep_eur_plink_pwas.sh
```

Use `PWAS_CHR_LIST` to restrict chromosomes:

```bash
PWAS_CHR_LIST="2 3 4" bash Analysis_code/pwas/4_1000Genome_keep_eur_plink_pwas.sh
```

## Step 5. Run S-MiXcan Association

Script:

```bash
Rscript Analysis_code/pwas/5_CARDMPRI_run_analysis_pwas.R
```

Inputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_input/chr*_mw_gwas_input_cardmpri_pwas.rds
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_pwas.bim
Results/pwas/cardmpri_workspace/cardmpri_filtered_id/filtered_chr*_hg38_012_pwas.raw
```

Outputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_chr*_result_pwas.csv
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_result_pwas.csv
```

For the real binary CARDMPRI analysis, provide case and control counts:

```bash
PWAS_N_CASES=<case_count> \
PWAS_N_CONTROLS=<control_count> \
Rscript Analysis_code/pwas/5_CARDMPRI_run_analysis_pwas.R
```

For technical smoke tests only, gaussian mode can be used without case/control
counts:

```bash
PWAS_CHR_LIST=1 PWAS_GWAS_FAMILY=gaussian \
Rscript Analysis_code/pwas/5_CARDMPRI_run_analysis_pwas.R
```

Step 5 result columns:

```text
gene_id
gene_name
chr
type
input_snp_num
Z_cardiomyocytes
p_cardiomyocytes
Z_other
p_other
p_join
```

## Step 6. Run Primo Pattern Inference

Script:

```bash
Rscript Analysis_code/pwas/6_run_Primo_pwas.R
```

Inputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_result_pwas.csv
Data/ensembl38.txt
```

Outputs:

```text
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_result_pwas_annotated.csv
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_table_pwas.csv
```

This step annotates genes with cytoband and runs:

```r
infer_celltype_patterns()
```

using cardiomyocyte and other p-values.

## Current Run Notes

As of the latest run:

```text
Training weights:
  Results/pwas/training_model_weights/weights_pwas_cardiomyocytes_other.csv
  17,615 weight rows
  6,048 unique genes
  4 skipped genes

Step 3 CARDMPRI merge:
  Results/pwas/cardmpri_workspace/cardmpri_input/
  Results/pwas/cardmpri_workspace/cardmpri_filtered_id/
  Complete for chromosomes 1-22
  Total merged rows: 16,197

Step 4:
  Tested successfully for chr1
  chr1 retained 1,296 variants and 633 EUR samples

Step 5:
  Tested successfully for chr1 in gaussian smoke-test mode
  chr1 result rows: 605
```

## Typical Full Workflow

```bash
# Step 2: train all weights, or train per chromosome and combine.
Rscript Analysis_code/pwas/2_train_model_pwas.R

# Step 3: merge weights with CARDMPRI GWAS.
Rscript Analysis_code/pwas/3_CARDMPRI_prepare_data_pwas.R

# Step 4: build LD reference files.
PWAS_PLINK2_ENV=plink2 PWAS_PLINK_ENV=plink2 \
bash Analysis_code/pwas/4_1000Genome_keep_eur_plink_pwas.sh

# Step 5: run CARDMPRI S-MiXcan association.
PWAS_N_CASES=<case_count> PWAS_N_CONTROLS=<control_count> \
Rscript Analysis_code/pwas/5_CARDMPRI_run_analysis_pwas.R

# Step 6: run Primo.
Rscript Analysis_code/pwas/6_run_Primo_pwas.R
```

