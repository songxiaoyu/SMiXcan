# safe_ACAT  calculates ACAT including catching error
#' @param p_values
#' @return ACAT p-value result
#' @expor
safe_ACAT <- function(p_values) {
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)  # Return NA or an alternative value
  })
}
# TO ADD: Z score beta, se_beta,... correlation updated
# change the name updated_MiXcan_association_test

# MiXcan_association_join
#' @param new_y gene expression level
#' @param new_cov covariates
#' @param new_outcome phenotype
#' @param family gaussian or binomial
#' @return p-value
#' @importFrom broom tidy
#' @importFrom dplyr filter mutate bind_cols
#' @importFrom magrittr %>%
#' @expor
MiXcan_association_join<-function (new_y, new_cov, new_outcome, family = gaussian)
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


# MiXcan_association_sep
#' @param new_y gene expression level
#' @param new_cov covariates
#' @param new_outcome phenotype
#' @param family gaussian or binomial
#' @return p-value
#' @importFrom broom tidy
#' @importFrom dplyr filter mutate bind_cols
#' @importFrom magrittr %>%
#' @expor
MiXcan_association_sep <- function(new_y, new_cov,
                                   new_outcome, family= gaussian){

  dat1 <- data.frame(cell_1 = new_y[,1], cov = new_cov,  y = new_outcome[,1])
  #
  dat2 <- data.frame(cell_2 = new_y[,2], cov = new_cov,  y = new_outcome[,1])
  if(identical(family, 'binomial')){
    dat1$y <- as.factor(dat1$y)
    dat2$y <- as.factor(dat2$y)
  }
  #
  # ft1 ft2
  ft1 <- glm(as.formula(paste0("y ~.")), data = dat1, family = family)
  ft2 <- glm(as.formula(paste0("y ~.")), data = dat2, family = family)

  res1_full <- broom::tidy(ft1)
  res1<- res1_full %>% filter(term == "cell_1")
  res2_full <- broom::tidy(ft2)
  res2 <- res2_full %>% filter(term == "cell_2")
  p_combined <- safe_ACAT(c(res1$p.value, res2$p.value))

  result <- res1 %>%
    dplyr::select(cell1_est = estimate,
                  cell1_se = std.error,
                  cell1_p = p.value) %>%
    bind_cols(res2 %>%
                dplyr::select(cell2_est = estimate,
                              cell2_se = std.error,
                              cell2_p = p.value)) %>%
    mutate(p_combined = p_combined) %>%
    as.data.frame()
  return(result)
}


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

