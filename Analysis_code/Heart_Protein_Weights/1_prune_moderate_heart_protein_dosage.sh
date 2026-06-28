#!/bin/bash
set -euo pipefail

# Create moderately LD-pruned GTEx EA dosage files for PWAS training.
#
# Motivation:
#   - Fully unpruned by_chr_nomiss has ~16M SNPs and is too large for local
#     training.
#   - The archived 500kb/r2=0.8 pruning is small but can make PWAS weights sparse.
#   - This script keeps more SNPs by default with 500kb/r2=0.95.
#
# Examples:
#   PWAS_CHR_LIST=1 Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh
#   PWAS_PRUNE_R2=0.9 Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh

PLINK="${PLINK:-/opt/anaconda3/envs/plink2/bin/plink2}"

PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
NOMISS_DIR="${PWAS_NOMISS_BY_CHR_DIR:-$PAPER_DIR/New generated files/codes/by_chr_nomiss}"

PRUNE_KB="${PWAS_PRUNE_KB:-500}"
PRUNE_STEP="${PWAS_PRUNE_STEP:-1}"
PRUNE_R2="${PWAS_PRUNE_R2:-0.95}"
TAG="${PRUNE_KB}kb_${PRUNE_STEP}_r2_${PRUNE_R2}"

OUTDIR="${PWAS_MODERATE_PRUNED_DIR:-$PAPER_DIR/New generated files/codes/moderate_pruned_by_chr_${TAG}}"
CHR_LIST="${PWAS_CHR_LIST:-$(printf "%s," {1..22})}"
CHR_LIST="${CHR_LIST%,}"

mkdir -p "$OUTDIR"

check_nonempty() {
  if [ ! -s "$1" ]; then
    echo "ERROR: missing or empty file: $1"
    exit 1
  fi
}

IFS=',' read -r -a CHRS <<< "$CHR_LIST"
for CHR in "${CHRS[@]}"; do
  echo "===== Moderate pruning chr${CHR}, ${TAG} ====="

  NOMISS_PREFIX="$NOMISS_DIR/GTEx_EA_chr${CHR}_nomiss"
  PRUNE_PREFIX="$OUTDIR/GTEx_EA_chr${CHR}_${TAG}"
  PRUNED_PREFIX="$OUTDIR/chr${CHR}_dosage_nomiss_LDpruned_${TAG}"
  PRUNED_RAW="${PRUNED_PREFIX}.raw"

  check_nonempty "${NOMISS_PREFIX}.pgen"
  check_nonempty "${NOMISS_PREFIX}.pvar"
  check_nonempty "${NOMISS_PREFIX}.psam"

  if [ -s "${PRUNE_PREFIX}.prune.in" ]; then
    echo "Prune list exists, skipping: ${PRUNE_PREFIX}.prune.in"
  else
    echo "Running LD pruning chr${CHR}"
    "$PLINK" \
      --pfile "$NOMISS_PREFIX" \
      --rm-dup exclude-all \
      --indep-pairwise "${PRUNE_KB}kb" "$PRUNE_STEP" "$PRUNE_R2" \
      --out "$PRUNE_PREFIX"
  fi
  check_nonempty "${PRUNE_PREFIX}.prune.in"

  if [ -s "$PRUNED_RAW" ]; then
    echo "Pruned raw exists, skipping: $PRUNED_RAW"
  else
    echo "Exporting moderate-pruned dosage chr${CHR}"
    "$PLINK" \
      --pfile "$NOMISS_PREFIX" \
      --extract "${PRUNE_PREFIX}.prune.in" \
      --export A \
      --out "$PRUNED_PREFIX"
  fi
  check_nonempty "$PRUNED_RAW"

  echo "===== Finished chr${CHR} ====="
done
