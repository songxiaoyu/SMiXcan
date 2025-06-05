
# regularized_inverse_cov calculates regularized inverse covariance
#' @param X X is the covariance matrix
#' @return a list of inverse matrix X_inv and regularized parameter lambda
#' @expor
regularized_inverse_cov <- function(X, lambda = NULL) {
  
  r <- abs(cor(X))
  if(r[1,2] > 0.9 || r[2,3] > 0.9){
    print('regularization applied')
    lambda = 0.3 * abs(r)^6
  }else{
    lambda = 0.3 * abs(r)^6
  }
  # Apply regularization if lambda > 0
  X_cor <- cov2cor(X)
  X_reg <- X_cor + lambda * diag(nrow(X))
  # X_cor_inv <- pseudoinverse(X_reg)
  X_cor_inv <- solve(X_reg)
  D <- diag(X)
  X_inv <- diag(1/sqrt(D)) %*% X_cor_inv %*% diag(1/sqrt(D))
  
  return(list(
    inv = X_inv,
    lambda = lambda
  ))
}


# run_S_MiXcan_r runs S-MiXcan with shrinkage (final version)
#' @param W1 weights for cell-type 1
#' @param W2 weights for cell-type 2
#' @param gwas_results gwas_results
#' @param cov_x_g covariance for genome
#' @param YtY Y is gene expression YtY = t(Y) %*% Y
#' @family0 binomial or gaussian
#' @Z0 Z-score for c(1,...,1)
#' @return a list of p-value from S-MiXcan  c(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
#' @expor
run_S_MiXcan_r <- function(W1, W2, gwas_results, cov_x_g, YtY, family0, Z_0=0){
  
  Beta <- gwas_results$Beta
  se_Beta <- gwas_results$se_Beta
  
  sig2_1g = unlist(t(matrix(W1)) %*% cov_x_g %*% matrix(W1))
  sig2_2g = unlist(t(matrix(W2)) %*% cov_x_g %*% matrix(W2))
  sig_l = sqrt(diag(cov_x_g))
  
  Z_1_sep = sum((W1* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_1g))
  Z_2_sep = sum((W2* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_2g))
  p_1_sep = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
  p_2_sep = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
  p_sep <- safe_ACAT(c(p_1_sep, p_2_sep))
  
  if(family0 == 'binomial'){
    Omega = diag(YtY)
    Omega0 = Omega[1]
    Omega1 = Omega[2]
    Omega2 = Omega[3]
    #cor = YtY[2,3]/(sqrt(Omega1)*sqrt(Omega2))
    corY = cor(YtY)
    if(abs(corY[1,2])>0.999999 || abs(corY[2,3]) >0.999999){
      print('cor==1')
      p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
    }else{
      S<-regularized_inverse_cov(YtY)$inv
      
      Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega0) * Z_0, sqrt(Omega1) *Z_1_sep, sqrt(Omega2)*Z_2_sep))
      p_0_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
      p_1_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_join[3]), lower.tail=FALSE)
    }
    
    p_join <- safe_ACAT(c(p_1_join, p_2_join))
    
    results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
  }else if(family0 == 'gaussian'){
    Omega = diag(YtY)
    Omega1 = Omega[1]
    Omega2 = Omega[2]
    
    cor = YtY[1,2]/(sqrt(Omega1)*sqrt(Omega2))
    if(is.na(cor)){
      p_1_join = NA
      p_2_join = NA
    } else if(abs(cor) == 1){
      print('cor==1')
      p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
    }else{
      rank_YtY <- qr(YtY)$rank
      dim_YtY <- nrow(YtY)
      if (rank_YtY < dim_YtY) {
        S <- ginv(YtY)
      } else {
        # Matrix is invertible
        S <- solve(YtY)
      }
      
      Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega1) *Z_1_sep, sqrt(Omega2)*Z_2_sep))
      p_1_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
    }
    
    p_join <- safe_ACAT(c(p_1_join, p_2_join))
    results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
    
  }else{
    print('family wrong')
  }
  return(results)
}


# run_S_MiXcan runs S-MiXcan without shrinkage 
#' @param W1 weights for cell-type 1
#' @param W2 weights for cell-type 2
#' @param gwas_results gwas_results
#' @param cov_x_g covariance for genome
#' @param YtY Y is gene expression YtY = t(Y) %*% Y
#' @family0 binomial or gaussian
#' @Z0 Z-score for c(1,...,1)
#' @return a list of p-value from S-MiXcan  c(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
#' @expor
run_S_MiXcan <- function(W1, W2, selected_snp, gwas_results, cov_x_g, YtY, family0, Z_0=0){
  
  Beta <- gwas_results$Beta
  se_Beta <- gwas_results$se_Beta
  
  sig2_1g = unlist(t(matrix(W1)) %*% cov_x_g %*% matrix(W1))
  sig2_2g = unlist(t(matrix(W2)) %*% cov_x_g %*% matrix(W2))
  sig_l = sqrt(diag(cov_x_g))
  
  Z_1_sep = sum((W1* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_1g))
  Z_2_sep = sum((W2* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_2g))
  
  p_1_sep = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
  p_2_sep = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
  p_sep <- safe_ACAT(c(p_1_sep, p_2_sep))
  
  if(family0 == 'binomial'){
    YtY=cov2cor(YtY)
    Omega = diag(YtY)
    Omega0 = Omega[1]
    Omega1 = Omega[2]
    Omega2 = Omega[3]
    cor = YtY[2,3]/(sqrt(Omega1)*sqrt(Omega2))
    if(abs(cor)==1){
      print('cor==1')
      p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
    }else{
      rank_YtY <- qr(YtY)$rank
      dim_YtY <- nrow(YtY)
      if (rank_YtY < dim_YtY) {
        S <- ginv(YtY)
      } else {
        # Matrix is invertible
        S <- solve(YtY)
      }
      Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega0)*Z_0, 
                                                          sqrt(Omega1)*Z_1_sep, 
                                                          sqrt(Omega2)*Z_2_sep))
      p_0_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
      p_1_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_join[3]), lower.tail=FALSE)
    }
    p_join <- safe_ACAT(c(p_1_join, p_2_join))
    
    results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
  }else if(family0 == 'gaussian'){
    Omega = diag(YtY)
    Omega1 = Omega[1]
    Omega2 = Omega[2]
    cor = YtY[1,2]/(sqrt(Omega1)*sqrt(Omega2))
    if(is.na(cor)){
      p_1_join = NA
      p_2_join = NA
    } else if(abs(cor) == 1){
      print('cor==1')
      p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
    }else{
      rank_YtY <- qr(YtY)$rank
      dim_YtY <- nrow(YtY)
        if (rank_YtY < dim_YtY) {
        S <- ginv(YtY)
      } else {
        # Matrix is invertible
        S <- solve(YtY)
      }
    
      Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega1) *Z_1_sep, sqrt(Omega2)*Z_2_sep))
      p_1_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
      p_2_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
    }
    
    p_join <- safe_ACAT(c(p_1_join, p_2_join))
    results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join)
    
  }else{
    print('family wrong')
  }
  return(results)
}


