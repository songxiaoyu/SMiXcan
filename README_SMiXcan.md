
# S-MiXcan

**S-MiXcan** is an R package for performing cell-type-aware **summary-based transcriptome-wide or proteome-wide association studies (TWAS/PWAS)**. It extends the MiXcan framework to analyze associations between genetically regulated gene expression (GReX) and traits using GWAS summary statistics, while accounting for cell-type-specific effects.

---

## 🔬 Overview

Traditional TWAS and PWAS models are based on bulk tissue expression data, potentially masking signals that are specific to particular cell types. **S-MiXcan** decomposes GReX into multiple cell-type-specific components and jointly tests their associations with phenotypic outcomes, enhancing biological resolution.

---

## 📦 Installation

You can install the development version of `S-MiXcan` from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install S-MiXcan from GitHub
devtools::install_github("YourUsername/SMiXcan")
```

> Replace `"YourUsername"` with your GitHub handle.

---

## 🧪 Example Usage

```r
library(SMiXcan)

# Load your trained elastic net models and cell-type proportion references
load("model_weights.RData")
load("cell_proportions.RData")

# Run S-MiXcan using GWAS summary statistics
results <- smixcan_assoc_test(
  gwas_data = "gwas_sumstats_chr1.txt",
  model_weights = model_weights,
  cell_props = cell_proportions,
  ld_reference = "1000G_EUR_chr1.ld",
  chr = 1
)

# View results
head(results)
```

---

## 📁 Folder Structure

```
SMiXcan/
├── R/                  # Core functions
├── data/               # Example models or data objects
├── man/                # Documentation
├── vignettes/          # Usage examples and tutorials
└── README.md           # Package description and instructions
```

---

## 📄 Citation

If you use **S-MiXcan** in your research, please cite:

> Zhu S, Song X, et al. (2025). *S-MiXcan: Cell-type-aware proteome-wide association analysis using GWAS summary statistics*. [Manuscript in preparation].

---

## 📫 Contact

For questions, please contact:

**Sinan Zhu**  
PhD Candidate, Duke-NUS Medical School  
Email: sinan.zhu@u.duke.nus.edu

---

## 🔒 License

This package is distributed under the MIT License. See `LICENSE` file for details.
