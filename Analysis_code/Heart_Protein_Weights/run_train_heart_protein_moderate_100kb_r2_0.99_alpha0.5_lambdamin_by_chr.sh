#!/bin/bash

# Train heart protein weights chromosome by chromosome with the selected less-sparse setup:
#   genotype input: moderate pruning, 100kb window, r2 = 0.99
#   MiXcan alpha: 0.5
#   glmnet lambda: lambda.min
#
# This writes one weight file per chromosome, then these can be combined with:
#   Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"
PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
GENO_DIR="${PWAS_GENO_RAW_DIR:-${PAPER_DIR}/New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99}"
LOG_DIR="${PWAS_MODERATE_TRAIN_LOG_DIR:-/private/tmp/pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin_logs}"
CHR_LIST="${PWAS_CHR_LIST:-$(printf "%s," {1..22})}"
CHR_LIST="${CHR_LIST%,}"
OUTPUT_PREFIX="${PWAS_OUTPUT_PREFIX:-_moderate_100kb_r2_0.99_alpha0.5_lambdamin}"

mkdir -p "${LOG_DIR}"

IFS=',' read -r -a CHRS <<< "${CHR_LIST}"

cd "${REPO_DIR}"

echo "Training heart protein weights with moderate pruning and lambda.min"
echo "Genotype dir: ${GENO_DIR}"
echo "Chromosomes: ${CHR_LIST}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Logs: ${LOG_DIR}"

for chr in "${CHRS[@]}"; do
  echo "===== Training chr${chr} ====="
  PWAS_GENO_RAW_DIR="${GENO_DIR}" \
  PWAS_CHR_FILTER="${chr}" \
  PWAS_OUTPUT_SUFFIX="${OUTPUT_PREFIX}_chr${chr}" \
  PWAS_MIXCAN_ALPHA=0.5 \
  PWAS_LAMBDA_CHOICE=min \
  PWAS_WEIGHT_EPS=1e-8 \
  "${RSCRIPT_BIN}" Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R \
    > "${LOG_DIR}/train_chr${chr}.log" 2>&1
  echo "===== Finished chr${chr} ====="
done

echo "All requested chromosomes finished."
echo "Combine outputs with:"
echo "  ${RSCRIPT_BIN} Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R"
