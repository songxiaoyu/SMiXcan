# SMiXcan

SMiXcan is an R package for cell-type-aware TWAS using trained SNP weights,
GWAS summary statistics, and an LD reference panel. The repository also includes
analysis workflows for the breast cell-type analyses and heart protein PWAS
analyses.

## Repository Layout

```text
R/                      R package source code
man/                    R package documentation
data/                   Small bundled example datasets
Analysis_code/          Reproducible analysis workflows
  2CellTypes/           Breast 2-cell-type workflow
  3CellTypes/           Breast 3-cell-type workflow
  Heart_Protein_Weights/ Heart protein weight training workflow
  HERMES/               HERMES DCM PWAS workflow
  HEARTFAIL/            HEARTFAIL PWAS workflow
```

Historical experiments, old PWAS scripts, generated results, and previous
nested package copies are not tracked in this repository.

## Installation

Install dependent packages:
```r
install.packages("devtools")
install_github("petraf01/BayesDeBulk/BayesDeBulk")

if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("limma")
install.packages("matrixStats")
install.packages("nnls")
install.packages("R.methodsS3")
install.packages("lcmix",repos="http://r-forge.r-project.org")
devtools::install_github("kjgleason/Primo")
```



Install our own package from GitHub :

```r
devtools::install_github("songxiaoyu/SMiXcan")
```





## Core Workflow

The package supports four main steps:

1. Estimate cell-type proportions from bulk expression.
2. Train cell-type-specific genetic prediction models.
3. Run SMiXcan association testing from GWAS summary statistics.
4. Optionally classify cell-type association patterns using PRIMO.

### 1. Estimate Cell-Type Fractions

```r
library(SMiXcan)

data(exprB_example, markers_example)

res_pi <- pi_estimation_K(
  exprB = exprB_example,
  markers = markers_example,
  seed = 1,
  n.iter = 1000,
  burn.in = 200
)

head(res_pi$cell_fraction)
```

### 2. Train Cell-Type-Specific Models

```r
data(x_example, y_example, pi_k)

set.seed(1)
foldid <- sample(rep(1:4, length.out = nrow(x_example)))

fit <- MiXcan_train_K(
  y = y_example,
  x = x_example,
  pi_k = pi_k,
  foldid = foldid,
  yName = "GENE_A1"
)

fit$type
fit$W
```

Key outputs include:

- `type`: model class, such as `CellTypeSpecific`, `NonSpecific`, or
  `NoPredictor`.
- `W`: SNP-by-cell-type weight matrix used for downstream association testing.
- `beta.SNP.by.cell`: SNP effects returned separately by cell type.
- `glmnet.cell` and `glmnet.tissue`: fitted penalized models.

### 3. Run SMiXcan Association Testing

```r
data(gwas_example)

res_assoc <- SMiXcan_assoc_test_K(
  W = fit$W,
  gwas_results = gwas_example,
  x_g = x_example,
  n0 = 1000,
  n1 = 1000,
  family = "binomial",
  regularization = "fixed",
  reg_scale = 0.1
)

res_assoc$p_join
res_assoc$p_join_vec
```

`SMiXcan_assoc_test_K()` filters SNPs with all-zero weights and supports two
regularization modes:

```r
regularization = "fixed"
regularization = "estimate"
```

Important outputs:

- `Z_join`, `p_join_vec`: regularized joint cell-type Z-scores and p-values.
- `p_join`: ACAT-combined p-value across cell types.
- `Z_sep`, `p_sep`: pre-regularization component statistics.
- `reg_scale_selected`: selected regularization scale.
- `reg_condition`: condition number after regularization.
- `mode`: `joint`, `separate`, or `empty`.

For most workflows, call the package function directly through:

```r
library(SMiXcan)
```

Do not source local copies of `R/S-MiXcan_K.R`.

### 4. PRIMO Pattern Classification

```r
data(merged_example)

res_primo <- infer_celltype_patterns(
  merged = merged_example,
  pvals_names = c("p_1_ct2", "p_2_ct2"),
  p_join_name = "p_join_ct2",
  type_col = "type_ct2",
  fdr_cutoff = 0.1
)

res_primo$tab_specific_patterns
```

For two cell types, PRIMO classifies non-null patterns such as:

- `10`: associated with cell type 1 only
- `01`: associated with cell type 2 only
- `11`: associated with both cell types

## Analysis Workflows

The `Analysis_code/` folder contains project-specific scripts. These scripts are
kept in GitHub for reproducibility but are not included in the R package build.

### Breast 2-Cell-Type Workflow

```text
Analysis_code/2CellTypes/
```

This is the original 2-cell-type breast workflow, including model training,
BCAC GWAS preparation, LD reference extraction, association testing, PRIMO, and
plotting.

### Breast 3-Cell-Type Workflow

```text
Analysis_code/3CellTypes/
```

This is the original 3-cell-type breast workflow.

### Heart Protein Weight Training

```text
Analysis_code/Heart_Protein_Weights/
```

This folder trains the heart protein prediction weights used by the HERMES and
HEARTFAIL PWAS workflows. The active weight model is documented in:

```text
Analysis_code/Heart_Protein_Weights/CURRENT_WEIGHTS_SOURCE.md
```

### HERMES PWAS Workflow

```text
Analysis_code/HERMES/
```

Main steps:

```bash
python Analysis_code/HERMES/1_liftover_hermes_gwas.py
Rscript Analysis_code/HERMES/2_HERMES_prepare_data_pwas.R
bash Analysis_code/HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh
Rscript Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
Rscript Analysis_code/HERMES/5_run_Primo_hermes_pwas.R
Rscript Analysis_code/HERMES/6_plot_hermes_pwas.R
```

HERMES default sample counts:

```text
HERMES_N_CASES=14256
HERMES_N_CONTROLS=1199156
```

Step 4 supports:

```bash
HERMES_ASSOC_REG_MODE=fixed
HERMES_ASSOC_REG_MODE=estimate
```

The active step 4 script intentionally keeps only `fixed` and `estimate`
regularization modes.

### HEARTFAIL PWAS Workflow

```text
Analysis_code/HEARTFAIL/
```

Main steps:

```bash
python Analysis_code/HEARTFAIL/1_liftover_heartfail_gwas.py
Rscript Analysis_code/HEARTFAIL/2_HEARTFAIL_prepare_data_pwas.R
bash Analysis_code/HEARTFAIL/3_1000Genome_keep_eur_plink_heartfail_pwas.sh
HEARTFAIL_N_CASES=<case_count> HEARTFAIL_N_CONTROLS=<control_count> \
  Rscript Analysis_code/HEARTFAIL/4_HEARTFAIL_run_analysis_pwas.R
Rscript Analysis_code/HEARTFAIL/5_run_Primo_heartfail_pwas.R
Rscript Analysis_code/HEARTFAIL/6_plot_heartfail_pwas.R
```

HEARTFAIL case/control counts must be set explicitly before step 4 unless the
script is edited for a specific GWAS release.

Step 4 supports:

```bash
HEARTFAIL_ASSOC_REG_MODE=fixed
HEARTFAIL_ASSOC_REG_MODE=estimate
```

## Notes For Development

- The repository now uses a single top-level R package. The previous nested
  `SMiXcan/` package copy is not part of the active repository layout.
- `Analysis_code/` is workflow code, not package source.
- The local source package can be rebuilt with `R CMD build .`.
- The current package association function should be called from
  `library(SMiXcan)` rather than by sourcing local R files.
