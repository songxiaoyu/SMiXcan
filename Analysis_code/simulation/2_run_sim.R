
safe_ACAT <- function(p_values) {
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)  # Return NA or an alternative value
  })
}

run_simulation <- function(Data_batch, sim_result_10, i){
  for (j in 1:B){
    sim_data = Data_batch[[j]]
    sim_data <- sim_center(sim_data)
    # p_results <- train_prediction_model(sim_data$y.train, sim_data$x.train,pi_example)
    p_results <- train_prediction_model(sim_data$y.train, sim_data$x.train, sim_data$pi.train)
    W1 <- p_results$W1
    W2 <- p_results$W2
    selected_snp <- p_results$selected_snp
    
    #run gwas
    gwas_results <- run_gwas(sim_data$x.test, sim_data$D.test, 'gaussian')
    
    # read covariance matrix and YTY
    x_g = sim_data$x.train[ ,selected_snp]
    cov_x_g <- cov(as.matrix(sim_data$x.train)) 
    
    XtX =t(sim_data$x.test) %*% as.matrix(sim_data$x.test)
    tilde_y_1 = sim_data$x.test %*% W1 
    tilde_y_2 = sim_data$x.test %*% W2 
    tilde_y_1 = tilde_y_1 - mean(tilde_y_1)
    tilde_y_2 = tilde_y_2 - mean(tilde_y_2)
    
    dat <- data.frame(tilde_y1 = tilde_y_1, tilde_y2 = tilde_y_2, tilde_D = sim_data$D.test[,1])
    MiXcan_prediction_result = data.frame(cell_1=tilde_y_1, cell_2=tilde_y_2)
    # MiXcan sep
    MiXcan_association_result_sep <- MiXcan_association_sep(new_y = MiXcan_prediction_result, 
                                                            new_cov = data.frame('cov'=rep(0,nrow(sim_data$x.test))), new_outcome = sim_data$D.test, family  = "gaussian")
    sim_result_10[(i-1)*B+j, 'p_m_sep_1'] = MiXcan_association_result_sep$cell1_p
    sim_result_10[(i-1)*B+j, 'p_m_sep_2'] = MiXcan_association_result_sep$cell2_p
    sim_result_10[(i-1)*B+j, 'p_m_sep'] = MiXcan_association_result_sep$p_combined
    
    # MiXcan join
    MiXcan_association_result_join <- MiXcan_association_join(new_y = MiXcan_prediction_result,
                                                                new_cov = data.frame('cov'=rep(0,nrow(sim_data$x.test))), new_outcome = sim_data$D.test, family  = "gaussian")
      
    sim_result_10[(i-1)*B+j, 'p_m_join_1'] = MiXcan_association_result_join$cell1_p
    sim_result_10[(i-1)*B+j, 'p_m_join_2'] = MiXcan_association_result_join$cell2_p
    sim_result_10[(i-1)*B+j, 'p_m_join'] = MiXcan_association_result_join$p_combined
    
    
    # S-MiXcan 
   
    YtY <- t(as.matrix(dat[, c('tilde_y1', 'tilde_y2')])) %*% as.matrix(dat[, c('tilde_y1', 'tilde_y2')])
    S_MiXcan_results <-run_S_MiXcan(W1, W2, selected_snp, gwas_results, cov_x_g, YtY, family0 = 'gaussian')

    sim_result_10[(i-1)*B+j, 'p_s_sep_1'] = S_MiXcan_results$p_1_sep
    sim_result_10[(i-1)*B+j, 'p_s_sep_2'] = S_MiXcan_results$p_2_sep
    sim_result_10[(i-1)*B+j, 'p_s_sep'] = S_MiXcan_results$p_sep
    sim_result_10[(i-1)*B+j, 'p_s_join_1'] = S_MiXcan_results$p_1_join
    sim_result_10[(i-1)*B+j, 'p_s_join_2'] = S_MiXcan_results$p_2_join
    sim_result_10[(i-1)*B+j, 'p_s_join'] = S_MiXcan_results$p_join
  }
  return(sim_result_10)
}
