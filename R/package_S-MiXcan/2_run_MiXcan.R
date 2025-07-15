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
