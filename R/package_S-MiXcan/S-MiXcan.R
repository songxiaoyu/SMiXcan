
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

# SMiXcan_assoc_test
# estimate YtY from reference genome
# run_S_MiXcan_r runs S-MiXcan with shrinkage (final version)
#' @param W1 a p by 1 matrix of weights for cell-type 1 (p is the number of snps in the neighborhood of the gene)
#' @param W2 a p by 1 matrix of weights for cell-type 2
#' @param gwas_results a list from gwas results analysis results: the first element is the vector of effect size Beta, and the second vector is standard error of Beta se_Beta
#' @param ref_g reference genome # name changed
#' @param ref_ld ld covariance matrix from reference # name changed
#' @n0 number of controls
#' @n1 number of cases
#' @family 'binomial' or 'gaussian'
#' @return a list of p-value from S-MiXcan  c(p_1_sep, p_2_sep, p_1_join, p_2_join, p_join)
#' @importFrom tibble lst
#' @export
SMiXcan_assoc_test <- function(W1, W2, gwas_results, x_g, ld_ref, n0, n1, family){
  y_0 <- rep(1,n1 + n0)
  D <- c(rep(1, n1), rep(0, n0))
  ft_0 <- glm(D ~ 0 + y_0, family = family)
  Z_0 <- coef(summary(ft_0))['y_0', 3]


  W = cbind(W1,W2)
  Y = x_g %*% W
  Y_scaled=scale(Y)
  Y <- cbind(rep(1,nrow(x_g)), Y_scaled)
  colnames(Y) = c('Y0', 'Y1', 'Y2')
  YtY = t(Y) %*% Y

  Beta <- gwas_results$Beta
  se_Beta <- gwas_results$se_Beta

  sig2_1g = unlist(t(matrix(W1)) %*% ld_ref %*% matrix(W1))
  sig2_2g = unlist(t(matrix(W2)) %*% ld_ref %*% matrix(W2))
  sig_l = sqrt(diag(ld_ref))

  Z_1_sep = sum((W1* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_1g))
  Z_2_sep = sum((W2* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_2g))
  p_1_sep = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
  p_2_sep = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
  p_sep <- safe_ACAT(c(p_1_sep, p_2_sep))

  Omega = diag(YtY)
  Omega0 = Omega[1]
  Omega1 = Omega[2]
  Omega2 = Omega[3]
  corY = cor(YtY)
  if(abs(corY[1,2])>0.999999 || abs(corY[2,3]) >0.999999){
    print('cor==1')
    p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
    Z_1_join = Z_1_sep
    p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
    Z_2_join = Z_2_sep
  }else{
    S<-regularized_inverse_cov(YtY)$inv
    Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega0) * Z_0, sqrt(Omega1) *Z_1_sep, sqrt(Omega2)*Z_2_sep))
    p_0_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
    p_1_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
    Z_1_join = Z_join[1]
    p_2_join = 2*pnorm(q=abs(Z_join[3]), lower.tail=FALSE)
    Z_2_join = Z_join[2]
  }

  p_join <- safe_ACAT(c(p_1_join, p_2_join))
  results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join, Z_1_join, Z_2_join)

  return(results)
}

run_S_MiXcan_r <- function(W1, W2, gwas_results, x_g, ld_ref, n0, n1, family){
  y_0 <- rep(1,n1 + n0)
  D <- c(rep(1, n1), rep(0, n0))
  ft_0 <- glm(D ~ 0 + y_0, family = family)
  Z_0 <- coef(summary(ft_0))['y_0', 3]

  W = cbind(W1,W2)
  Y = x_g %*% W
  Y_scaled=scale(Y)
  Y <- cbind(rep(1,nrow(x_g)), Y_scaled)
  colnames(Y) = c('Y0', 'Y1', 'Y2')
  YtY = t(Y) %*% Y

  Beta <- gwas_results$Beta
  se_Beta <- gwas_results$se_Beta

  sig2_1g = unlist(t(matrix(W1)) %*% ld_ref %*% matrix(W1))
  sig2_2g = unlist(t(matrix(W2)) %*% ld_ref %*% matrix(W2))
  sig_l = sqrt(diag(ld_ref))

  Z_1_sep = sum((W1* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_1g))
  Z_2_sep = sum((W2* Beta) * sig_l / se_Beta) * (1 / sqrt(sig2_2g))
  p_1_sep = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
  p_2_sep = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
  p_sep <- safe_ACAT(c(p_1_sep, p_2_sep))

  Omega = diag(YtY)
  Omega0 = Omega[1]
  Omega1 = Omega[2]
  Omega2 = Omega[3]
  corY = cor(YtY)
  if(abs(corY[1,2])>0.999999 || abs(corY[2,3]) >0.999999){
      print('cor==1')
      p_1_join = 2*pnorm(q=abs(Z_1_sep), lower.tail=FALSE)
      Z_1_join = Z_1_sep
      p_2_join = 2*pnorm(q=abs(Z_2_sep), lower.tail=FALSE)
      Z_2_join = Z_2_sep
    }else{
      S<-regularized_inverse_cov(YtY)$inv
      Z_join <-diag(1/sqrt((diag(S)))) %*% S %*% matrix(c(sqrt(Omega0) * Z_0, sqrt(Omega1) *Z_1_sep, sqrt(Omega2)*Z_2_sep))
      p_0_join = 2*pnorm(q=abs(Z_join[1]), lower.tail=FALSE)
      p_1_join = 2*pnorm(q=abs(Z_join[2]), lower.tail=FALSE)
      Z_1_join = Z_join[1]
      p_2_join = 2*pnorm(q=abs(Z_join[3]), lower.tail=FALSE)
      Z_2_join = Z_join[2]
    }

    p_join <- safe_ACAT(c(p_1_join, p_2_join))
    results<- lst(p_1_sep, p_2_sep, p_sep, p_1_join, p_2_join, p_join, Z_1_join, Z_2_join)

  return(results)
}
