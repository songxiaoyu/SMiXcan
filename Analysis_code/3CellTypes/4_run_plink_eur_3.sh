#!/bin/bash

# Explicitly source conda setup
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate plink2

# Script name retained for history. The kept 3pi workflow below uses the
# current pi3 filenames already present in the shared workspace.

# Base paths
CHR_LIST=$(echo {1..22})
DATA_DIR="/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene"
GWAS_ID_DIR="/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/bcac2020_filtered_id"

for chr in $CHR_LIST; do
  conda activate plink2
  echo "Processing chr${chr}..."

  # Decompress files if needed (only once per chr)
  # Added -f to force overwrite if file exists to prevent errors on re-runs
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
      plink2 --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
      plink2 --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  # Extract the SNPs used by the kept 3pi workflow.
  plink2 \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/bcac2020_filtered_chr${chr}_gwas_id_pi3.txt" \
    --make-bed \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_pi3"

  # Convert to 012 format (additive coding)
  conda activate plink
  plink \
    --bfile "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_pi3" \
    --recodeA \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_012_pi3"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."
