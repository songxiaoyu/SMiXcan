# HEARTFAIL PWAS Analysis

This folder contains the HEARTFAIL GWAS version of the PWAS/S-MiXcan workflow.

## Input

Default GWAS file:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes.tsv.bgz
```

The raw GWAS columns are expected to include:

```text
variant, minor_allele, beta, se, pval
```

`variant` is parsed as `chr:pos:ref:alt` in hg19/GRCh37. The script uses
`minor_allele` as the GWAS effect allele, so if `minor_allele == REF`, the
baseline allele is `ALT`; if `minor_allele == ALT`, the baseline allele is
`REF`.

## Output

Default workspace:

```text
Results/heartfail_pwas/heartfail_workspace
```

Main outputs:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_table_pwas.csv
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ.pdf
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ_split_venn.pdf
```

## Workflow

### Step 1. Lift HEARTFAIL GWAS from hg19 to hg38

```bash
python3 Analysis_code/HEARTFAIL/1_liftover_heartfail_gwas.py
```

This writes:

```text
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes_hg38_rsid_chr_pos.txt
```

### Step 2. Merge PWAS weights with HEARTFAIL GWAS

```bash
Rscript Analysis_code/HEARTFAIL/2_HEARTFAIL_prepare_data_pwas.R
```

Optional overrides:

```bash
HEARTFAIL_CHR_LIST=1
HEARTFAIL_GWAS_HG38_FILE=/path/to/lifted.tsv.gz
PWAS_WEIGHTS_FILE=/path/to/weights.csv
HEARTFAIL_WORKSPACE_DIR=/path/to/workspace
```

### Step 3. Build EUR LD reference files

```bash
PLINK2_BIN=/opt/anaconda3/envs/plink2/bin/plink2 \
PLINK_BIN=/opt/anaconda3/envs/plink/bin/plink \
bash Analysis_code/HEARTFAIL/3_1000Genome_keep_eur_plink_heartfail_pwas.sh
```

Use the correct local PLINK/PLINK2 paths for your machine.

### Step 4. Run S-MiXcan association

You must provide HEARTFAIL case/control counts before running the default
binomial model:

```bash
HEARTFAIL_N_CASES=<case_count> \
HEARTFAIL_N_CONTROLS=<control_count> \
Rscript Analysis_code/HEARTFAIL/4_HEARTFAIL_run_analysis_pwas.R
```

Regularization options:

```bash
HEARTFAIL_ASSOC_REG_MODE=fixed
HEARTFAIL_ASSOC_REG_SCALE=0.005
```

or:

```bash
HEARTFAIL_ASSOC_REG_MODE=estimate
```

Estimate-mode internal defaults are handled by the `SMiXcan` package function.

### Step 5. Run Primo pattern inference

```bash
Rscript Analysis_code/HEARTFAIL/5_run_Primo_heartfail_pwas.R
```

### Step 6. Plot QQ and split Venn

```bash
Rscript Analysis_code/HEARTFAIL/6_plot_heartfail_pwas.R
```

## Notes

The HEARTFAIL case/control counts are not hard-coded because the GWAS file does
not encode them unambiguously. Confirm the phenotype counts before step 4.
