
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
  if (is.na(res1$p.value* res2$p.value)) {
    cat("error")
    p_combined=NA
  }else{
    p_combined <- ACAT::ACAT(c(res1$p.value, res2$p.value))
  }
  
  
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
