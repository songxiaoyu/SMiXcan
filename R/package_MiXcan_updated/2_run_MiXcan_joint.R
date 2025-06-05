
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
  if (cor(new_y[, 1], new_y[, 2]) == 1) {
    res2 <- res1
  }
  if (cor(new_y[, 1], new_y[, 2]) == -1) {
    res2 <- res1 %>% mutate(estimate = -estimate)
  }
  if(is.na(res1$p.value) | is.na(res2$p.value)){
    p_combined <- NA
  }else{
  p_combined <- ACAT::ACAT(c(res1$p.value, res2$p.value))
  }
  result <- res1 %>% dplyr::select(cell1_est = estimate, cell1_se = std.error, 
                                   cell1_p = p.value) %>% bind_cols(res2 %>% dplyr::select(cell2_est = estimate, 
                                                                                           cell2_se = std.error, cell2_p = p.value)) %>% mutate(p_combined = p_combined) %>% 
    as.data.frame()
  return(result)
}


