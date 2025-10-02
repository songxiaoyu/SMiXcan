import os
import pandas as pd
from pyliftover import LiftOver

# Initialize LiftOver once
lo = LiftOver('hg19', 'hg38')

# LiftOver function
def convert_position(row):
    result = lo.convert_coordinate(row['chr.iCOGs'], row['Position.iCOGs'])
    return result[0][1] if result else None

# Directory containing input files
input_dir = "/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/BCAC/Breast_Cancer_Risk_2020"
output_dir = input_dir  # Or use a subdirectory like f"{input_dir}/lifted"

# Chromosomes to process
# chroms = list(range(1, 23)) + ['X']  # chr1 to chr22 + chrX

chroms = list(range(10, 12))
for chrom in chroms:
    input_file = os.path.join(input_dir, f"chr{chrom}_icogs.txt")
    output_file = os.path.join(output_dir, f"chr{chrom}_icogs_hg38.csv")

    print(f"Processing {input_file}...")

    try:
        # Load data
        df = pd.read_csv(input_file, delim_whitespace=True, header=0, low_memory=False, on_bad_lines='skip')

        # Standardize chr format
        df['chr.iCOGs'] = df['chr.iCOGs'].apply(lambda x: f'chr{x}')
        df['Position.iCOGs'] = pd.to_numeric(df['Position.iCOGs'], errors='coerce')

        # Apply liftover
        df['POS_hg38'] = df.apply(convert_position, axis=1)
        df = df.dropna(subset=['POS_hg38'])

        # Save result
        df.to_csv(output_file, index=False, na_rep="")

        print(f"✓ Saved {len(df)} lifted positions to {output_file}\n")

    except Exception as e:
        print(f"⚠️ Error processing chr{chrom}: {e}\n")
