# train_prediction_model
#' @param y.train gene expression level in training data 
#' @param x.train genome in training data 
#' @param pi.train pi estimate in training data 
#' @return a list of W1 W2 and selected_snp
#' @expor
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
