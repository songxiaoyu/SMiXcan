#' Compute a regularized inverse covariance matrix
#'
#' This helper stabilizes covariance-matrix inversion by applying a simple
#' correlation-based shrinkage step before inversion.
#'
#' @name regularized_inverse_cov
#' @title Regularized Covariance Matrix Inversion
#'
#' @param X Numeric covariance matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{inv}{The regularized inverse matrix.}
#'   \item{lambda}{The regularization parameter used.}
#' }
#'
#' @export

regularized_inverse_cov <- function(X) {

  r <- abs(stats::cor(X))
  #lambda = 0.001 * abs(r)^6
  lambda = 0.1 * abs(r)^2
  #lambda = 0.01 * abs(r)^2
  #lambda = 0.005 * abs(r)^3
  #lambda = 0.1 * abs(r)^3
  # Apply regularization if lambda> 0
  X_cor <- stats::cov2cor(X)
  X_reg <- X_cor + lambda * diag(nrow(X))
  X_cor_inv <- solve(X_reg)
  D <- diag(X)
  X_inv <- diag(1/sqrt(D)) %*% X_cor_inv %*% diag(1/sqrt(D))

  return(list(
    inv = X_inv,
    lambda = lambda
  ))
}

#' Run the K-cell-type S-MiXcan association test
#'
#' This function combines cell-type-specific prediction weights, GWAS summary
#' statistics, and a reference genotype panel to compute per-cell-type joint
#' association statistics.
#' @name SMiXcan_assoc_test_K
#' @title  S-MiXcan Association Test in K cell types
#'
#' @param W Numeric p-by-K matrix of cell-type-specific SNP weights, where
#'   p is the number of SNPs and K is the number of cell types.
#' @param gwas_results List containing GWAS summary statistics with elements
#'   \code{Beta} and \code{se_Beta}.
#' @param x_g Reference-panel genotype matrix with individuals in rows and
#'   SNPs in columns.
#' @param n0 Number of controls.
#' @param n1 Number of cases.
#' @param family Either \code{"binomial"} or \code{"gaussian"}.
#' @return A list containing:
#' \describe{
#'   \item{Z_join}{Vector of per-cell-type joint Z-scores.}
#'   \item{p_join_vec}{Vector of per-cell-type joint p-values.}
#'   \item{p_join}{Combined ACAT p-value across the K cell types.}
#' }
#'
#' @importFrom tibble lst
#' @export
SMiXcan_assoc_test_K <- function(W,
                                 gwas_results,
                                 x_g,
                                 n0,
                                 n1,
                                 family = c("binomial", "gaussian")) {
  family <- match.arg(family)

  # ---- Input checks ----
  Beta    <- as.numeric(gwas_results$Beta)
  se_Beta <- as.numeric(gwas_results$se_Beta)

  if (!is.matrix(W))
    stop("W must be a numeric matrix of dimension p * K.")

  p <- nrow(W)
  K <- ncol(W)

  if (ncol(x_g) != p)
    stop("Number of SNPs (columns) in x_g must match nrow(W).")

  if (length(Beta) != p || length(se_Beta) != p)
    stop("Length of Beta and se_Beta must match nrow(W).")

  # ---- LD and SNP variances ----
  cov_x <- stats::var(x_g)           # p * p LD (covariance) matrix
  sig_l <- sqrt(diag(cov_x))         # per-SNP SD

  # ---- Separate component Z-scores for each cell type ----
  Z_sep <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    wk <- W[, k]
    sig2_gk <- drop(t(wk) %*% cov_x %*% wk)  # Var(predicted expression)
    if (sig2_gk > 0) {
      num <- sum((wk * Beta) * sig_l / se_Beta)
      Z_sep[k] <- num / sqrt(sig2_gk)
    } else {
      Z_sep[k] <- NA_real_
    }
  }

  p_sep <- ifelse(
    is.na(Z_sep),
    NA_real_,
    2 * stats::pnorm(abs(Z_sep), lower.tail = FALSE)
  )

  # Start from the separate-model result and overwrite if the joint step is stable.
  Z_join     <- Z_sep
  p_join_vec <- p_sep
  mode       <- "separate"

  # ---- Joint test ----
  if (family == "binomial") {
    # 1. Fit the null intercept-only logistic model.
    D   <- c(rep(1L, n1), rep(0L, n0))
    fit0 <- stats::glm(D ~ 1, family = stats::binomial())
    coef0 <- summary(fit0)$coefficients
    Z0    <- coef0["(Intercept)", "z value"]

    # 2. Build predicted expression for each cell type.
    Yhat <- x_g %*% W             # n * K
    Y_scaled <- scale(Yhat)       # standardize columns
    Y <- cbind(1, Y_scaled)       # [Intercept, cell-type 1..K]
    colnames(Y) <- c("Y0", paste0("Y", seq_len(K)))

    YtY   <- crossprod(Y)         # (K+1) * (K+1)
    Omega <- diag(YtY)

    # Correlation matrix on the predictor scale.
    corY <- stats::cov2cor(YtY)

    # Skip the joint step when the predictor correlation is too close to singular.
    if (!any(is.na(corY)) &&
        !any(abs(corY[upper.tri(corY)]) > 0.999999)) {

      S <- regularized_inverse_cov(YtY)$inv
      v <- c(
        sqrt(Omega[1]) * Z0,
        sqrt(Omega[-1]) * Z_sep
      )

      Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
      Z_join <- as.numeric(Z_full[-1])  # drop intercept
      p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
      mode <- "joint"
    }

  } else if (family == "gaussian") {
    # Gaussian joint step uses only the cell-type predictors.
    Yhat <- x_g %*% W
    Y_scaled <- scale(Yhat)             # n * K
    colnames(Y_scaled) <- paste0("Y", seq_len(K))

    YtY   <- crossprod(Y_scaled)        # K * K
    Omega <- diag(YtY)
    corY  <- stats::cov2cor(YtY)

    if (!any(is.na(corY)) &&
        !any(abs(corY[upper.tri(corY)]) > 0.999999)) {

      S <- regularized_inverse_cov(YtY)$inv
      v <- sqrt(Omega) * Z_sep          # length K

      Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
      Z_join <- as.numeric(Z_full)
      p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
      mode <- "joint"
    }

  } else {
    stop("family must be 'binomial' or 'gaussian'")
  }

  # ---- ACAT-combined p-value across K cell types ----
  p_join <- safe_ACAT(p_join_vec)

  list(
    Z_join    = Z_join,
    p_join_vec = p_join_vec,
    p_join    = p_join
  )
}

