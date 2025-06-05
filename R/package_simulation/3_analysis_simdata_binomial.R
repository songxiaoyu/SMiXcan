library(data.table)
library(lme4)
# library('secr')
library(glmnet)
library(doParallel)
library(doRNG)
library(ACAT)
#library("VariantAnnotation")
#library('fourPNO')
library(lme4)
library(MiXcan)
library(tibble) 
library(tidyr)
library(dplyr)
library(MASS)



setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder')
file_names_1 <- list.files(path = "code_MiXcan_updated", pattern = "\\.R$")
setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_MiXcan_updated')
lapply(file_names_1, source)

setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder')
file_names_2 <- list.files(path = "code_S-MiXcan", pattern = "\\.R$")
setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_S-MiXcan')
lapply(file_names_2, source)

setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder')
file_names_3 <- list.files(path = "code_Simulation", pattern = "\\.R$")
#file_names_3 <- file_names_3[2]
setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_Simulation')
lapply(file_names_3[2], source)




X_pool_filtered <- read.csv('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_Simulation/X_pool_filtered.csv')


B = 100 # batch 
ITR = 2 # total iterations
eta <- matrix(c(0, 0, 0.2, 0.2, 0, 0.2, 0.2, 0, -0.2, 0.2), ncol = 2, byrow = TRUE)

for (nonzero_beta2 in c(2)){
for(t in 1:5){
  sim_result_bin = data.frame(matrix(ncol = 12, nrow = ITR*B))
  colnames(sim_result_bin) <- c('p_s_sep_1','p_s_sep_2','p_s_sep','p_m_sep_1','p_m_sep_2','p_m_sep',
                               'p_s_join_1','p_s_join_2','p_s_join','p_m_join_1','p_m_join_2','p_m_join')
  for(i in 1:ITR){
    Data_batch = DataGen_newhap_binom(mc=B, n.train=300, n.test=3000, p=50, b0=0, 
                                      nonzero_beta1=1, nonzero_beta2=nonzero_beta2 , gammas=eta[t, ], 
                                      var1=1, var2=1, seed=(i+1),group = 'heter xy', X_pool=X_pool_filtered,snp_region = 1:50)
    filename = paste0('Data/Data_batch_0313_heter_binom_nonzero_beta2_',nonzero_beta2,'_', t,'setting',i,'.RData')
    saveRDS(Data_batch,filename)
     Data_batch = readRDS(filename)
    print(i)
    for(j in 1:B){
      sim_data = Data_batch[[j]]
      sim_data <- sim_center_binom(sim_data)
      
      # train prediction model
      p_results <- train_prediction_model(sim_data$y.train, sim_data$x.train, sim_data$pi.train)
      W1 <- p_results$W1
      W2 <- p_results$W2
      selected_snp <- p_results$selected_snp
      
      #run gwas
      gwas_results <- run_gwas(sim_data$x.test, sim_data$D.test, 'binomial')
      
      # read covariance matrix and YTY
      x_g = sim_data$x.train[ ,selected_snp]
      cov_x_g <- cov(as.matrix(sim_data$x.train)) 
      
      tilde_y_1 = sim_data$x.test %*% W1 
      tilde_y_2 = sim_data$x.test %*% W2 
      tilde_y_0 <- rep(1,nrow(sim_data$x.test))
      dat <- data.frame(tilde_y0 =tilde_y_0, tilde_y1 = tilde_y_1 - mean(tilde_y_1), tilde_y2 = tilde_y_2 - mean(tilde_y_2), tilde_D = sim_data$D.test[,1])
      YtY <- t(as.matrix(dat[, c('tilde_y0', 'tilde_y1', 'tilde_y2')])) %*% as.matrix(dat[, c('tilde_y0', 'tilde_y1', 'tilde_y2')])
      
      # MiXcan sep
      MiXcan_prediction_result = data.frame(cell_1=tilde_y_1- mean(tilde_y_1), cell_2=tilde_y_2- mean(tilde_y_2))
      # MiXcan_prediction_result = data.frame(cell_1=tilde_y_1, cell_2=tilde_y_2)
      
      MiXcan_association_result_sep <- MiXcan_association_sep(new_y = MiXcan_prediction_result, 
                                                              new_cov = data.frame('cov'=rep(0,nrow(sim_data$x.test))), new_outcome = sim_data$D.test, family  = "binomial")
      sim_result_bin[(i-1)*B+j, 'p_m_sep_1'] = MiXcan_association_result_sep$cell1_p
      sim_result_bin[(i-1)*B+j, 'p_m_sep_2'] = MiXcan_association_result_sep$cell2_p
      sim_result_bin[(i-1)*B+j, 'p_m_sep'] = MiXcan_association_result_sep$p_combined
      
      # MiXcan join
      if(is.na(cor(tilde_y_1, tilde_y_2))| cor(tilde_y_1, tilde_y_2) == 1){
        sim_result_bin[(i-1)*B+j, 'p_m_join'] = NA
        next
      }
      else{
        MiXcan_association_result_join <- MiXcan_association_join(new_y = MiXcan_prediction_result,
                                                                  new_cov = data.frame('cov'=rep(0,nrow(sim_data$x.test))), new_outcome = sim_data$D.test, family  = 'binomial')
        
        sim_result_bin[(i-1)*B+j, 'p_m_join_1'] = MiXcan_association_result_join$cell1_p
        sim_result_bin[(i-1)*B+j, 'p_m_join_2'] = MiXcan_association_result_join$cell2_p
        sim_result_bin[(i-1)*B+j, 'p_m_join'] = MiXcan_association_result_join$p_combined
        
        # S-MiXcan 
        ft_0 <- glm(tilde_D ~ 0 + tilde_y0, data = dat, family = 'binomial')
        Z_0 <- coef(summary(ft_0))['tilde_y0', 3]
        S_MiXcan_results <-run_S_MiXcan(W1, W2, selected_snp, gwas_results, cov_x_g, YtY, family0 = 'binomial', Z_0=Z_0)
        
        p_1_sep = S_MiXcan_results$p_1_sep
        p_2_sep = S_MiXcan_results$p_2_sep
        
        if(is.na(p_1_sep) | is.na(p_2_sep)){
          sim_result_bin[(i-1)*B+j, ] = NA
          next
        }
        sim_result_bin[(i-1)*B+j, 'p_s_sep_1'] = S_MiXcan_results$p_1_sep
        sim_result_bin[(i-1)*B+j, 'p_s_sep_2'] = S_MiXcan_results$p_2_sep
        sim_result_bin[(i-1)*B+j, 'p_s_sep'] = S_MiXcan_results$p_sep
        
        sim_result_bin[(i-1)*B+j, 'p_s_join_1'] = S_MiXcan_results$p_1_join
        sim_result_bin[(i-1)*B+j, 'p_s_join_2'] = S_MiXcan_results$p_2_join
        sim_result_bin[(i-1)*B+j, 'p_s_join'] = S_MiXcan_results$p_join
        print((i-1)*B+j)
      }
    }
  }
  filename = paste0('Data/result_batch_0313_heter_binom_nonzero_beta2_',nonzero_beta2,'_', t,'setting',i,'.csv')
  write.csv(sim_result_bin, filename)
}
}
