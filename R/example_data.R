#' Example GWAS, LD Covariance, and Genotype Reference Data (chr21)
#'
#' Example inputs for S-MiXcan analysis on chromosome 21.
#'
#' @format
#' \describe{
#'   \item{filtered_mw_gwas_input}{A data frame of GWAS summary statistics
#'         merged with cell-type-specific prediction weights.}
#'   \item{cov_ref}{A numeric symmetric matrix giving LD covariance
#'         between SNPs (rows and columns are SNP rsIDs).}
#'   \item{X_ref}{A numeric matrix of reference genotypes (rows = samples,
#'         columns = SNP rsIDs).}
#' }
#'
#' @examples
#' data(example_data_chr21)
#' head(filtered_mw_gwas_input)
#' dim(cov_ref)
#' dim(X_ref)
"filtered_mw_gwas_input"
"cov_ref"
"X_ref"
