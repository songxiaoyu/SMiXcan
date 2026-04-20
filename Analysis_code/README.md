# Analysis Pipeline Guide

This folder contains the analysis scripts used for the `SMiXcan` paper workflows. The code is organized into three parts:

- `2CellTypes/`: two-cell-type analysis in breast tissue
- `3CellTypes/`: three-cell-type analysis in breast tissue

The scripts were originally written for a local analysis environment with several large reference datasets stored outside GitHub. This README explains the intended run order, what each script does, and which files need to exist before each stage can run.

## Before You Start

There are two kinds of inputs in this project:

- Small-to-moderate files that are now expected in Dropbox under `Paper_SMiXcan/Data` or `Paper_SMiXcan/Results`
- Very large local resources that are still referenced from absolute paths outside Dropbox, such as GTEx expression/genotype files, BCAC GWAS files, and 1000 Genomes reference genotype files

In other words, this repository is the code layer, but not a full self-contained data bundle.

## Expected Folder Layout

The current scripts assume the following Dropbox layout:

```text
Paper_SMiXcan/
├── Data/
│   ├── BCAC/
│   │   ├── chr1_icogs_hg38.csv
│   │   ├── ...
│   │   └── chr22_icogs_hg38.csv
│   ├── DRIVE/
│   ├── GTEx/
│   └── ensembl38.txt
├── Results/
│   ├── BayesDeBulk_pi_3ct_GTEx.tsv
│   ├── SMiXcanK_results/
│   └── drive_result_full_lam_new.csv
├── Figure/
└── Github/
    └── Analysis_code/
```

The code also expects several large local working directories that are not mirrored into Dropbox:

- `.../code_RealData/RealData/GTEx_Data`
- `.../2pi`
- `.../3pi`
- local 1000 Genomes PLINK reference files

If you run the scripts on a different machine, the first thing to check is the hard-coded path section near the top of each script.

## Package Setup

Install the package from the repository root before running the analysis scripts:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "ggplot2",
  "ggrepel", "cowplot", "ggforce", "bacon", "glmnet", "lme4",
  "doParallel", "doRNG", "janitor", "tibble", "DBI", "RSQLite",
  "MASS", "ACAT", "Primo"
))

install.packages(".", repos = NULL, type = "source")
```

Many scripts use `library(SMiXcan)`, so the package should be installed before running them.

## Two-Cell-Type Pipeline

The two-cell analysis lives in [2CellTypes](./2CellTypes).

### Overview

This workflow:

1. trains a two-cell-type prediction model in GTEx breast tissue
2. evaluates agreement between MiXcan and S-MiXcan in the DRIVE dataset
3. prepares BCAC summary statistics for chromosome-wise analysis
4. builds European-only LD reference subsets from 1000 Genomes
5. runs the S-MiXcan association scan across all chromosomes
6. annotates results and performs Primo pattern inference
7. generates the main figures

### Script Order

1. [2CellTypes/1_pi_estimate_2.R](./2CellTypes/1_pi_estimate_2.R)

This file is only a note. The two-cell proportions used here are derived in the training workflow from the 3-cell BayesDeBulk estimates, so there is no standalone 2-cell `pi` estimation script to run.

2. [2CellTypes/2_train_model_2.R](./2CellTypes/2_train_model_2.R)

Purpose:
- loads GTEx breast expression, covariates, genotype, and PredictDB elastic-net SNP sets
- loads `Results/BayesDeBulk_pi_3ct_GTEx.tsv`
- collapses the 3-cell proportions into epithelial versus other
- trains the two-cell MiXcan model gene by gene

Main output:
- `Paper_SMiXcan/Results/weights_miXcan_full_pi2.csv`

3. [2CellTypes/3_DRIVE_Analysis_gtex.R](./2CellTypes/3_DRIVE_Analysis_gtex.R)

Purpose:
- loads the trained two-cell weights
- loads DRIVE genotype and phenotype data
- runs single-SNP GWAS in DRIVE
- compares individual-level MiXcan and summary-statistics S-MiXcan

Main output:
- `Paper_SMiXcan/Results/drive_result_full_lam_new.csv`

Notes:
- this script contains a small local helper called `smixcan_assoc_for_drive()`
- that helper is only for this comparison script and is not part of the package API

4. [2CellTypes/4_BCAC_prepare_data.R](./2CellTypes/4_BCAC_prepare_data.R)

Purpose:
- reads chromosome-wise BCAC GWAS summary statistics
- matches BCAC SNPs to the trained two-cell weight table
- writes chromosome-specific merged inputs

Input location:
- `Paper_SMiXcan/Data/BCAC/chr*_icogs_hg38.csv`

Source:
- BCAC GWAS summary statistics were downloaded from the BCAC summary results page: `https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations`

Main outputs:
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_filtered_id/`
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_input/`

5. [2CellTypes/5_1000Genome_keep_eur_plink.sh](./2CellTypes/5_1000Genome_keep_eur_plink.sh)

Purpose:
- uses the SNP lists from step 4
- extracts European-only reference genotypes from 1000 Genomes with PLINK

Main output location:
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_filtered_id/`

This step requires PLINK and local 1000 Genomes reference files.

6. [2CellTypes/6_BCAC2020_run_analysis_2.R](./2CellTypes/6_BCAC2020_run_analysis_2.R)

Purpose:
- loads the chromosome-wise merged BCAC inputs
- loads the European-only LD reference panels from step 5
- runs `SMiXcan_assoc_test_K()` gene by gene across chromosomes
- combines all chromosome outputs

Main output location:
- `Paper_SMiXcan/Results/bcac2020_result/`

7. [2CellTypes/7_run_Primo.R](./2CellTypes/7_run_Primo.R)

Purpose:
- reads `bcac2020_result_pi2.csv`
- adds gene name and cytoband annotation using `Data/ensembl38.txt`
- runs `infer_celltype_patterns()`
- prepares the supplementary table for the two-cell analysis

Main outputs in `Paper_SMiXcan/Results/SMiXcanK_results`:
- `bcac2020_result_pi2_annotated.csv`
- `tableS2.csv`

8. [2CellTypes/8_plot.R](./2CellTypes/8_plot.R)

Purpose:
- uses the annotated two-cell BCAC output and the DRIVE comparison output
- generates the main manuscript figures

Main outputs in `Paper_SMiXcan/Figure`:
- `Figure1_ABC_scatter.pdf`
- `Figure2_BC.pdf`
- `Final_Figure_ABC.pdf`

## Three-Cell-Type Pipeline

The three-cell analysis lives in [3CellTypes](./3CellTypes).

### Overview

This workflow:

1. estimates three-cell-type proportions in GTEx breast tissue
2. trains a three-cell MiXcan model
3. prepares BCAC chromosome-wise inputs
4. builds European-only LD reference subsets
5. runs the three-cell S-MiXcan scan
6. annotates the results and runs Primo
7. generates the supplementary figures

### Script Order

1. [3CellTypes/1_pi_estimation_gtex_3.R](./3CellTypes/1_pi_estimation_gtex_3.R)

Purpose:
- estimates cell fractions with `pi_estimation_K()`
- uses a marker list for adipose/endothelial, fibroblast, and epithelial compartments

Main output:
- `Paper_SMiXcan/Results/BayesDeBulk_pi_3ct_GTEx.tsv`

2. [3CellTypes/2_train_model_3.R](./3CellTypes/2_train_model_3.R)

Purpose:
- loads GTEx expression, covariates, genotype, and PredictDB elastic-net SNP sets
- loads the three-cell BayesDeBulk proportions
- trains a three-cell MiXcan model gene by gene

Main output:
- `Paper_SMiXcan/Results/weights_miXcan_full_pi3.csv`

3. [3CellTypes/3_BCAC2020_prepare_data_3.R](./3CellTypes/3_BCAC2020_prepare_data_3.R)

Purpose:
- reads chromosome-wise BCAC GWAS summary statistics
- matches BCAC SNPs to the three-cell weight table
- writes chromosome-specific merged inputs

Input location:
- `Paper_SMiXcan/Data/BCAC/chr*_icogs_hg38.csv`

Source:
- BCAC GWAS summary statistics were downloaded from the BCAC summary results page: `https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations`

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/`

Note:
- these archived Dropbox workspace folders contain the available 3-cell chromosome-specific merged inputs and SNP ID files
- the original workspace folder keeps the name `baca2020_input/` in Dropbox because that is how it was named in the local run directory

4. [3CellTypes/4_run_plink_eur_3.sh](./3CellTypes/4_run_plink_eur_3.sh)

Purpose:
- uses the SNP lists from step 3
- extracts the 1000 Genomes reference genotypes used by the kept `3pi` workflow with PLINK

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/bcac2020_filtered_id/`

5. [3CellTypes/5_BCAC2020_run_analysis_3.R](./3CellTypes/5_BCAC2020_run_analysis_3.R)

Purpose:
- loads the chromosome-wise merged BCAC inputs
- loads the European-only LD reference panels from step 4
- runs `SMiXcan_assoc_test_K()` gene by gene across chromosomes
- combines all chromosome outputs

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/bcac2020_result/`

6. [3CellTypes/6_BCACresult.R](./3CellTypes/6_BCACresult.R)

Purpose:
- reads the combined three-cell result table
- adds gene and cytoband annotation
- runs `infer_celltype_patterns()`
- prepares the supplementary table for the three-cell analysis

Main outputs in `Paper_SMiXcan/Results/SMiXcanK_results`:
- `bcac2020_result_pi3_annotated.csv`
- `tableS3.csv`

7. [3CellTypes/7_plotS1.R](./3CellTypes/7_plotS1.R)

Purpose:
- uses the annotated three-cell BCAC output and the DRIVE comparison output
- generates the supplementary figures

Main outputs in `Paper_SMiXcan/Figure`:
- `FigureS1_111.pdf`
- `FigureS2_111.pdf`

## Typical Run Order for Reproducing the Main Results

If you want the main paper outputs rather than every intermediate step, the practical run order is:

1. install the `SMiXcan` package from this repository
2. make sure Dropbox `Data`, `Results`, and `Figure` folders exist
3. make sure the large local GTEx, BCAC, and 1000 Genomes resources are available
4. run `2CellTypes/2_train_model_2.R`
5. run `2CellTypes/3_DRIVE_Analysis_gtex.R`
6. run `2CellTypes/4_BCAC_prepare_data.R`
7. run `2CellTypes/5_1000Genome_keep_eur_plink.sh`
8. run `2CellTypes/6_BCAC2020_run_analysis_2.R`
9. run `2CellTypes/7_run_Primo.R`
10. run `2CellTypes/8_plot.R`

For the three-cell supplementary analysis, run the corresponding `3CellTypes` scripts in numeric order.

## Important Path Caveats

There are still a few path assumptions to be aware of:

- Several scripts use absolute local paths under `/Users/.../Downloads/S-MiXcan_code_folder/...`
- The two-cell pipeline currently points to a local `2pi` workspace for large intermediate files
- The three-cell pipeline currently points to the kept `3pi` workspace
- The archived 3-cell workspace retains one historical folder typo, `baca2020_input/`, and the scripts follow that name so they match the stored run artifacts

If your local machine uses a different folder name than `3pi`, update the related 3-cell scripts consistently before running.

## Notes for New Users

- The scripts are intended to be run one file at a time, not as a single automated pipeline
- Most intermediate files are chromosome-specific and are written into the local `2pi` or 3-cell working directories
- Final cleaned outputs are written into Dropbox `Results`
- Final figures are written into Dropbox `Figure`
- The file [Document for downloading reference genome.docx](./Document%20for%20downloading%20reference%20genome.docx) contains the download/setup notes for the external reference resources used in this project, including the GTEx-related local data setup and the reference genome resources used by the PLINK steps

