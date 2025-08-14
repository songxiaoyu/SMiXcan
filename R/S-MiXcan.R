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
#'   \item{Z_1_join}{Z-score for cell-type 1 (joint model)}
#'   \item{p_1_join}{P-value for cell-type 1 (joint model)}
#'   \item{Z_2_join}{Z-score for cell-type 2 (joint model)}
#'   \item{p_2_join}{P-value for cell-type 2 (joint model)}
#'   \item{p_join}{Combined p-value from joint model (ACAT)}
#' }
#'
#' @importFrom tibble lst
#' @export

SMiXcan_assoc_test <-  function (W1, W2, gwas_results, x_g, ld_ref, n0, n1, family = c("binomial","gaussian"))
{
  family <- match.arg(family)

  # ---- Inputs ----
  Beta    <- as.numeric(gwas_results$Beta)
  se_Beta <- as.numeric(gwas_results$se_Beta)

  # Basic length checks (silent alignment assumptions will bite later)
  stopifnot(length(W1) == length(W2),
            length(W1) == ncol(ld_ref),
            length(Beta) == length(W1),
            length(se_Beta) == length(W1))

  # ---- Separate-component Z ----
  sig2_1g <- drop(t(W1) %*% ld_ref %*% W1)
  sig2_2g <- drop(t(W2) %*% ld_ref %*% W2)
  sig_l   <- sqrt(diag(ld_ref))

  Z_1_sep <- if (sig2_1g > 0) sum((W1 * Beta) * sig_l / se_Beta) / sqrt(sig2_1g) else NA_real_
  Z_2_sep <- if (sig2_2g > 0) sum((W2 * Beta) * sig_l / se_Beta) / sqrt(sig2_2g) else NA_real_

  p_1_sep <- if (is.na(Z_1_sep)) NA_real_ else 2 * pnorm(abs(Z_1_sep), lower.tail = FALSE)
  p_2_sep <- if (is.na(Z_2_sep)) NA_real_ else 2 * pnorm(abs(Z_2_sep), lower.tail = FALSE)
  p_sep   <- safe_ACAT(c(p_1_sep, p_2_sep))

  # ---- Joint test(s) ----
  p_1_join <- p_2_join <- p_join <- NA_real_
  Z_1_join <- Z_2_join <- NA_real_

  if (family == "binomial") {
    # Intercept-only logistic model to get Z0
    D   <- c(rep(1L, n1), rep(0L, n0))
    fit <- glm(D ~ 1, family = binomial())
    Z_0 <- coef(summary(fit))["(Intercept)", "z value"]

    W <- cbind(W1, W2)
    Yhat <- x_g %*% W
    Y_scaled <- scale(Yhat)                          # standardize columns
    Y <- cbind(1, Y_scaled)                          # [Intercept, Y1, Y2]
    colnames(Y) <- c("Y0", "Y1", "Y2")

    YtY   <- crossprod(Y)
    Omega <- diag(YtY)
    # Correlation among columns of Y (use cov2cor of YtY to be consistent)
    corY <- cov2cor(YtY)

    if (any(is.na(corY)) ||
        abs(corY[1, 2]) > 0.999999 || abs(corY[1, 3]) > 0.999999 || abs(corY[2, 3]) > 0.999999) {
      # Fall back to separate if near-singular
      p_1_join <- p_1_sep; Z_1_join <- Z_1_sep
      p_2_join <- p_2_sep; Z_2_join <- Z_2_sep
    } else {
      S <- regularized_inverse_cov(YtY)$inv  # user helper
      Z_join <- diag(1 / sqrt(diag(S))) %*% S %*%
        matrix(c(sqrt(Omega[1]) * Z_0,
                 sqrt(Omega[2]) * Z_1_sep,
                 sqrt(Omega[3]) * Z_2_sep), ncol = 1)

      # Indices: 1=intercept, 2=component 1, 3=component 2
      Z_1_join <- as.numeric(Z_join[2])
      Z_2_join <- as.numeric(Z_join[3])
      p_1_join <- 2 * pnorm(abs(Z_1_join), lower.tail = FALSE)
      p_2_join <- 2 * pnorm(abs(Z_2_join), lower.tail = FALSE)
    }

  } else if (family == "gaussian") {
    # No intercept in the joint step (two components)
    W <- cbind(W1, W2)
    Yhat <- x_g %*% W
    Y_scaled <- scale(Yhat)
    colnames(Y_scaled) <- c("Y1","Y2")

    YtY   <- crossprod(Y_scaled)
    Omega <- diag(YtY)
    corY  <- cov2cor(YtY)

    if (any(is.na(corY)) || abs(corY[1, 2]) > 0.999999) {
      p_1_join <- p_1_sep
      p_2_join <- p_2_sep
      Z_1_join <- Z_1_sep
      Z_2_join <- Z_2_sep
    } else {
      S <- regularized_inverse_cov(YtY)$inv
      Z_join <- diag(1 / sqrt(diag(S))) %*% S %*%
        matrix(c(sqrt(Omega[1]) * Z_1_sep,
                 sqrt(Omega[2]) * Z_2_sep), ncol = 1)

      Z_1_join <- as.numeric(Z_join[1])
      Z_2_join <- as.numeric(Z_join[2])
      p_1_join <- 2 * pnorm(abs(Z_1_join), lower.tail = FALSE)
      p_2_join <- 2 * pnorm(abs(Z_2_join), lower.tail = FALSE)
    }
  } else {
    stop("family must be 'binomial' or 'gaussian'")
  }

  p_join <- safe_ACAT(c(p_1_join, p_2_join))

  # Return a base R list with names your caller expects
  list(
    Z_1_join = Z_1_join,
    p_1_join = p_1_join,
    Z_2_join = Z_2_join,
    p_2_join = p_2_join,
    p_join   = p_join
  )
}

