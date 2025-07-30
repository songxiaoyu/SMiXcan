#' Regularized Inverse of a Covariance Matrix
#'
#' Computes a regularized inverse of a covariance matrix using correlation shrinkage.
#'
#' @name regularized_inverse_cov
#' @title Regularized Covariance Matrix Inversion
#'
#' @param X A numeric covariance matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{inv}{The regularized inverse matrix.}
#'   \item{lambda}{The regularization parameter used.}
#' }
#'
#' @export

regularized_inverse_cov <- function(X) {

  r <- abs(cor(X))
  lambda = 0.3 * abs(r)^6
  # Apply regularization if lambda > 0
  X_cor <- cov2cor(X)
  X_reg <- X_cor + lambda * diag(nrow(X))
  X_cor_inv <- solve(X_reg)
  D <- diag(X)
  X_inv <- diag(1/sqrt(D)) %*% X_cor_inv %*% diag(1/sqrt(D))

  return(list(
    inv = X_inv,
    lambda = lambda
  ))
}

#' S-MiXcan Association Test with Shrinkage
#'
#' Performs the S-MiXcan association test using GWAS summary statistics and reference LD data,
#' applying shrinkage-based regularization to stabilize inverse covariance estimation.
#'
#' @param W1 A p-by-1 matrix of weights for cell-type 1 (where p is the number of SNPs in the gene region).
#' @param W2 A p-by-1 matrix of weights for cell-type 2.
#' @param gwas_results A list containing GWAS summary statistics, with components \code{Beta} and \code{se_Beta}.
#' @param x_g A genotype matrix for the reference panel (individuals × SNPs).
#' @param ld_ref A p-by-p LD covariance matrix from the reference panel.
#' @param n0 Number of controls.
#' @param n1 Number of cases.
#' @param family Either \code{"binomial"} or \code{"gaussian"} (used for fitting the null model).
#'
#' @return A list containing:
#' \describe{
#'   \item{p_1_sep}{P-value for cell-type 1 (separate model)}
#'   \item{p_2_sep}{P-value for cell-type 2 (separate model)}
#'   \item{p_sep}{Combined p-value from separate models (ACAT)}
#'   \item{p_1_join}{P-value for cell-type 1 (joint model)}
#'   \item{p_2_join}{P-value for cell-type 2 (joint model)}
#'   \item{p_join}{Combined p-value from joint model (ACAT)}
#'   \item{Z_1_join}{Z-score for cell-type 1 (joint model)}
#'   \item{Z_2_join}{Z-score for cell-type 2 (joint model)}
#' }
#'
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

