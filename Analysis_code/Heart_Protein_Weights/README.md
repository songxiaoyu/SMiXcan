# Heart Protein Weights Training

This folder contains only the scripts used to train GTEx heart protein
prediction weights for downstream HERMES and HEARTFAIL PWAS/S-MiXcan analyses.

The outputs are saved under:

`Results/heart_protein_weights/training_model_weights`

## Current Weights

The current retained model is:

`weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv`

It uses two components:

1. Cardiomyocytes
2. Other cell types combined

Training settings:

- Genotype dosage: GTEx EA heart dosage after moderate pruning
- Pruning: 100 kb window, `r2 = 0.99`
- MiXcan alpha: `0.5`
- glmnet lambda: `lambda.min`

See `CURRENT_WEIGHTS_SOURCE.md` for the full source record.

## Scripts

1. `1_prune_moderate_heart_protein_dosage.sh`
   - Creates moderately pruned dosage files from no-missing GTEx dosage.

2. `2_train_heart_protein_weights.R`
   - Trains per-chromosome protein prediction weights.

3. `run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh`
   - Runs step 2 across chromosomes with the selected moderate-pruning setup.

4. `3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R`
   - Combines per-chromosome weights into one full weights table.

## Reproduce Current Weights

```bash
bash Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh
bash Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh
Rscript Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R
```

Downstream HERMES and HEARTFAIL prepare scripts now default to:

`Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv`

Set `PWAS_WEIGHTS_FILE` only when intentionally testing a different weights file.
