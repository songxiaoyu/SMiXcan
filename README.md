
# S-MiXcan

**S-MiXcan** is an R package for performing cell-type-aware **summary-based transcriptome-wide association studies (TWAS)**. It extends the MiXcan framework to analyze associations between genetically regulated gene expression (GReX) and traits using GWAS summary statistics, while accounting for cell-type-specific effects.

---

## 🔬 Overview

Traditional transcriptome-wide association studies (TWAS) predict gene expression at the tissue level and test associations of predicted tissue-level expression with disease, often disregarding the diverse roles of different cell types in disease etiology. MiXcan is a recent approach that enables cell-type-aware TWAS but requires individual-level genotype data. In this study, we introduce S-MiXcan, a summary statistics-based, cell-type-aware TWAS framework. S-MiXcan utilizes summary statistics to infer associations between diseases and predicted expression at the cell-type level. Simulation studies demonstrate that S-MiXcan maintains well-controlled type I error rates and yields results comparable to those obtained using individual genotype data.

---

## 📦 Installation

You can install the development version of `S-MiXcan` from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install S-MiXcan from GitHub
devtools::install_github("songxiaoyu/SMiXcan")
```

---

## 🧪 Example Usage

```r
library(SMiXcan)

# Load example data
data(example_data_chr21)  
# Provides: filtered_mw_gwas_input, cov_ref, X_ref

# Split data by gene
split_df <- split(filtered_mw_gwas_input, filtered_mw_gwas_input$gene_name)
filtered_list <- list()
for (gene in names(split_df)) {
  gene_df <- split_df[[gene]]
  W1 <- gene_df$weight_cell_1
  W2 <- gene_df$weight_cell_2
  W  <- cbind(W1, W2)
  selected_snp <- gene_df[W1 != 0 | W2 != 0, ]
  filtered_list[[gene]] <- list(W = W, selected_snp = selected_snp)
}

# Case/control counts 
n_case <- 133384
n_control <- 113789 + 18908

# Run S-MiXcan for each gene
G <- length(filtered_list)
real_result <- data.frame(matrix(ncol = 8, nrow = G))
colnames(real_result) <- c("gene_name", "gene_id", "type",
                           "Z_s_join_1","p_s_join_1",
                           "Z_s_join_2","p_s_join_2","p_s_join")

for (g in seq_len(G)) {
  gene <- names(split_df)[g]
  cat("Processing gene:", gene, "\n")

  W <- filtered_list[[gene]]$W
  selected_snp_id <- filtered_list[[gene]]$selected_snp$rsid

  gwas_results <- list(
    Beta    = filtered_list[[gene]]$selected_snp$beta.Gwas,
    se_Beta = filtered_list[[gene]]$selected_snp$SE.Gwas
  )

  cov_x_g <- cov_ref[selected_snp_id, selected_snp_id, drop = FALSE]
  x_g     <- X_ref[, selected_snp_id, drop = FALSE]

  # Perform association test
  S_MiXcan_results <- SMiXcan_assoc_test(
    W[,1], W[,2],
    gwas_results,
    x_g, cov_x_g,
    n_control, n_case,
    family = "binomial"
  )

  real_result[g, ] <- c(
    gene,
    filtered_mw_gwas_input[1, "gene_id"],
    filtered_mw_gwas_input[1, "type"],
    S_MiXcan_results$Z_1_join,
    S_MiXcan_results$p_1_join,
    S_MiXcan_results$Z_2_join,
    S_MiXcan_results$p_2_join,
    S_MiXcan_results$p_join
  )
}

# Inspect results
head(real_result)
```

---

## 📁 Folder Structure

```
SMiXcan/
├── R/                  # Core functions
├── data/               # Example data objects
├── man/                # Documentation
├── vignettes/          # Usage examples and tutorials
└── README.md           # Package description and instructions
```

---

## 📄 Citation

If you use **S-MiXcan** in your research, please cite this page

---

## 📫 Contact

For questions, please contact:

**Sinan Zhu**  
PhD Candidate, Duke-NUS Medical School  
Email: sinan.zhu@u.duke.nus.edu

---

## 🔒 License


