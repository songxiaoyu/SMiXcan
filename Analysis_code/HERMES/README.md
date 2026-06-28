# HERMES DCM PWAS Analysis

This folder contains scripts for running the PWAS/S-MiXcan workflow on the HERMES
DCM GWAS summary statistics.

## Input

Local HERMES folder:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Heart/HERMES
```

GWAS summary table:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR.tsv.gz
```

GWAS build:

```text
hg19 / GRCh37
```

Sample counts:

```text
cases    = 14,256
controls = 1,199,156
```

Relevant GWAS columns:

```text
chr, pos_b37, A1, A2, A1_beta, se, pval
```

`A1_beta` is interpreted as the effect size for `A1`. Downstream scripts use
`A1` as `Effect.Gwas` and `A2` as `Baseline.Gwas`.

## Outputs

Main output directory:

```text
Results/hermes_pwas/hermes_workspace
```

Figures:

```text
Figure/hermes_pwas
```

## Workflow

Run from the repository root.

### 1. Liftover HERMES GWAS from hg19 to hg38

```bash
python3 Analysis_code/HERMES/1_liftover_hermes_gwas.py
```

Expected output:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR_hg38_rsid.tsv.gz
```

### 2. Prepare HERMES PWAS input

```bash
Rscript Analysis_code/HERMES/2_HERMES_prepare_data_pwas.R
```

This creates:

```text
Results/hermes_pwas/hermes_workspace/hermes_input
Results/hermes_pwas/hermes_workspace/hermes_filtered_id
```

### 3. Extract 1000Genome EUR LD reference SNPs

Set PLINK paths/environments if needed, then run:

```bash
PWAS_PLINK2_ENV=plink2 \
PLINK2_BIN=plink2 \
PLINK_BIN=plink \
bash Analysis_code/HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh
```

### 4. Run HERMES PWAS/S-MiXcan association

```bash
Rscript Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
```

Defaults:

```text
HERMES_N_CASES=14256
HERMES_N_CONTROLS=1199156
HERMES_GWAS_FAMILY=binomial
HERMES_ASSOC_REG_MODE=fixed
HERMES_ASSOC_REG_SCALE=0.1
```

`HERMES_ASSOC_REG_MODE` supports `fixed` and `estimate`. Estimate-mode internal
defaults are handled by the `SMiXcan` package function.

### 5. Run Primo pattern annotation

```bash
Rscript Analysis_code/HERMES/5_run_Primo_hermes_pwas.R
```

### 6. Plot HERMES PWAS results

```bash
Rscript Analysis_code/HERMES/6_plot_hermes_pwas.R
```

## Notes

- This workflow uses the existing PWAS training weights by default:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other.csv
```

- To use a new retrained PWAS weight table, set:

```bash
PWAS_WEIGHTS_FILE=/path/to/new_weights.csv
```

- To run a subset of chromosomes for testing, set:

```bash
HERMES_CHR_LIST=1
```
