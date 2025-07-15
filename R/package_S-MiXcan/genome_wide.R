
# run_gwas calculates gwas results
#' @param X X is the genome matrix
#' @param D is the phenotype
#' @family0 binomial or gaussian
#' @return a list of inverse matrix X_inv and regularized parameter lambda
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
