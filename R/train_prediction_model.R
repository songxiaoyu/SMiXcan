#' Train Cell-Type-Specific Prediction Model
#'
#' Trains a cell-type-aware prediction model using MiXcan, estimating weights for SNPs across cell types.
#'
#' @name train_prediction_model
#' @title Train Cell-Type-Specific Prediction Model Using MiXcan
#'
#' @param y.train Gene expression levels in training data.
#' @param x.train Genotype matrix in training data.
#' @param pi.train Estimated cell-type proportions in training data.
#'
#' @return A list containing:
#' \describe{
#'   \item{W1}{Weight vector for cell type 1}
#'   \item{W2}{Weight vector for cell type 2}
#'   \item{selected_snp}{Indices of selected SNPs}
#' }
#'
#' @importFrom MiXcan MiXcan MiXcan_extract_weight
#' @export

train_prediction_model <- function(y.train, x.train, pi.train){
  foldid_eg <- sample(1:10, length(y.train), replace=T)
  n_train = nrow(x.train)
  MiXcan_result <- MiXcan(y.train, x.train,
                          cov = data.frame('cov'=rep(0,n_train)),
                          pi.train,
                          foldid = foldid_eg)

  # To get training weights
  W1 = matrix(MiXcan_result$beta.SNP.cell1[,2],ncol=1)
  W2 = matrix(MiXcan_result$beta.SNP.cell2[,2],ncol=1)

  MiXcan_weight_result <- MiXcan_extract_weight(model = MiXcan_result)
  selected_snp = as.numeric(substr(MiXcan_weight_result$xNameMatrix, 4, nchar(MiXcan_weight_result$xNameMatrix)))
  results <- lst(W1, W2, selected_snp)
  return(results)
}
