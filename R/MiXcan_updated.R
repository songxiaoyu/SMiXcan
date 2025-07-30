#' Safe Wrapper for ACAT P-value Combination
#'
#' This function computes a combined p-value using the ACAT method while handling possible numerical or input errors.
#' If an error occurs during computation, the function returns \code{NA} and prints an informative message.
#'
#' @param p_values A numeric vector of p-values to be combined using the ACAT method.
#'
#' @return A single numeric value: the combined p-value (or \code{NA} if ACAT fails).
#'
#' @importFrom ACAT ACAT
#' @export
safe_ACAT <- function(p_values) {
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)
  })
}

#' Cell-Type-Aware Association (Joint Model)
#'
#' This function fits a generalized linear model (GLM) for joint association testing of two cell-type-specific predicted expressions with phenotype,
#' adjusting for covariates. It uses ACAT to combine p-values from both predictors.
#'
#' @param new_y A matrix or data frame containing predicted expression levels for two cell types (columns named "cell_1" and "cell_2").
#' @param new_cov A data frame of covariates (e.g., age, principal components).
#' @param new_outcome A vector or one-column matrix of phenotype values (e.g., binary case/control or continuous trait).
#' @param family A GLM family object or string: either \code{"gaussian"} for continuous outcomes or \code{"binomial"} for binary outcomes.
#'
#' @return A data frame containing estimates, standard errors, and p-values for both cell types, as well as the combined p-value using ACAT.
#'
#' @importFrom broom tidy
#' @importFrom dplyr filter mutate select bind_cols
#' @importFrom stats as.formula cor glm
#' @importFrom magrittr %>%
#' @export
MiXcan_assoc_test<-function (new_y, new_cov, new_outcome, family = 'gaussian')
{
  dat <- data.frame(new_y, new_cov, y = new_outcome[,1])
  if(identical(family, 'binomial')){
    dat$y <- as.factor(dat$y)
  }
  ft <- glm(as.formula(paste0("y ~.")), data = dat, family = family)
  res <- broom::tidy(ft)
  res1 <- res %>% filter(term == "cell_1")
  res2 <- res %>% filter(term == "cell_2")

  cory <- tryCatch({
    cor(new_y[, 1], new_y[, 2])
  }, error = function(e) {
    message("Correlation failed: ", e$message)
    return(NA)
  })
  if (is.na(cory)){
    p_combined <- NA
  }else if (cory == 1) {
    res2 <- res1
  }else if (cory == -1) {
    res2 <- res1 %>% mutate(estimate = -estimate)
  }

  p_combined <- safe_ACAT(c(res1$p.value, res2$p.value))

  result <- res1 %>% dplyr::select(cell1_est = estimate, cell1_se = std.error,
                                   cell1_p = p.value) %>% bind_cols(res2 %>% dplyr::select(cell2_est = estimate,
                                                                                           cell2_se = std.error, cell2_p = p.value)) %>% mutate(p_combined = p_combined) %>%
    as.data.frame()
  return(result)
}


#' @title Cell-Type-Aware Association with Shrinkage (Ridge Logistic Regression)
#'
#' @description
#' This function performs cell-type-aware association analysis using ridge-regularized logistic regression
#' followed by inference using an unpenalized model. It estimates the individual and combined p-values for
#' predicted gene expressions (e.g., from cell type 1 and cell type 2).
#'
#' @param new_outcome A binary phenotype vector (0/1 or factor).
#' @param new_y A matrix or data frame of predicted expression for two cell types (columns = cell_1 and cell_2).
#' @param new_cov A data frame of covariates (optional). Can include categorical variables (e.g., country, study).
#' @param lambda Optional regularization parameter. If NULL, estimated based on correlation structure.
#'
#' @return A data frame containing effect size estimates, standard errors, and p-values for both cell types,
#'         and a combined p-value from ACAT.
#'
#' @importFrom glmnet glmnet
#' @importFrom stats model.matrix glm binomial na.omit
#' @importFrom broom tidy
#' @importFrom dplyr select bind_cols mutate
#' @importFrom magrittr %>%
#' @export

MiXcan_association_shrink <- function(new_outcome, new_y, new_cov = NULL, lambda = NULL) {
  # Combine everything into one dataframe
  df <- data.frame(
    outcome = new_outcome,
    new_y,
    new_cov,
    check.names = TRUE
  )

  # Remove rows with any missing values
  df <- na.omit(df)

  # Build model matrix for glmnet
  X_glmnet <- model.matrix(~ . - outcome -1, data = df)
  y_clean <- df$outcome

  # Compute lambda if needed
  if (is.null(lambda)) {
    if (ncol(X_glmnet) > 1) {
      r <- cor(X_glmnet)
      lambda <- 0.3 * max(abs(r[lower.tri(r)]))^6
      if (is.na(lambda)) lambda <- 1e-6
    } else {
      lambda <- 1e-6
    }
  }

  # Fit ridge logistic regression
  fit_ridge <- glmnet(X_glmnet, y_clean, family = "binomial", alpha = 0, lambda = lambda, standardize = TRUE)

  # Get nonzero coefficients (excluding intercept)
  coef_nonzero <- coef(fit_ridge)
  idx_nonzero <- which(coef_nonzero != 0)[-1]  # drop intercept

  if (length(idx_nonzero) < 1) {
    return(data.frame(
      cell1_est = NA,
      cell1_se = NA,
      cell1_p = NA,
      cell2_est = NA,
      cell2_se = NA,
      cell2_p = NA,
      p_combined = NA
    ))
  }

  # Instead of reusing X_glmnet column names,
  # always rebuild formula with the ORIGINAL variable names:
  var_names <- c(colnames(new_y), colnames(new_cov))

  formula_str <- paste(
    "outcome ~",
    paste(sprintf("`%s`", var_names), collapse = " + ")
  )
  formula_glm <- as.formula(formula_str)

  # Refit unpenalized glm on the original data
  fit_glm <- glm(formula_glm, data = df, family = binomial())

  # Extract summary
  s <- summary(fit_glm)$coefficients

  # Names of the first two predictors
  cell1_name <- colnames(new_y)[1]
  cell2_name <- colnames(new_y)[2]

  # Helper
  get_coef_info <- function(var_name) {
    if (var_name %in% rownames(s)) {
      c(
        estimate = s[var_name, "Estimate"],
        se = s[var_name, "Std. Error"],
        p = s[var_name, "Pr(>|z|)"]
      )
    } else {
      c(estimate=NA, se=NA, p=NA)
    }
  }

  res1 <- get_coef_info(cell1_name)
  res2 <- get_coef_info(cell2_name)

  # Combine p-values with Fisher's method
  p_values <- c(res1["p"], res2["p"])
  p_values <- p_values[!is.na(p_values)]

  if (length(p_values) == 0) {
    p_combined <- NA
  } else {
    p_combined <- safe_ACAT(p_values)}

  # Return result
  result <- data.frame(
    cell1_est = res1["estimate"],
    cell1_se = res1["se"],
    cell1_p = res1["p"],
    cell2_est = res2["estimate"],
    cell2_se = res2["se"],
    cell2_p = res2["p"],
    p_combined = p_combined
  )

  return(result)
}

