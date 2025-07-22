#' Run GWAS for Each SNP
#'
#' Performs univariate GWAS by regressing the phenotype on each SNP separately.
#'
#' @name run_gwas
#' @title Run Genome-Wide Association Study
#'
#' @param X A genotype matrix (individuals × SNPs).
#' @param D A phenotype vector (numeric or binary).
#' @param family0 A GLM family, either "binomial" or "gaussian".
#'
#' @return A list with:
#' \describe{
#'   \item{Beta}{Vector of effect sizes for each SNP}
#'   \item{se_Beta}{Vector of standard errors for each SNP}
#' }
#'
#' @importFrom tibble lst
#' @export

run_gwas <- function(X, D, family0){
  Beta = rep(0, ncol(X))
  se_Beta = rep(0, ncol(X))
  for(j in 1: ncol(X)){
    gwas_result<- glm(D ~ X[,j], family=family0)
    Beta[j] = gwas_result$coefficients[2]
    se_Beta[j] = coef(summary(gwas_result))[, "Std. Error"][2]
  }
  results <- lst(Beta, se_Beta)
  return(results)
}
