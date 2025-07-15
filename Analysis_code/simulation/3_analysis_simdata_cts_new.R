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
lapply(file_names_3[3], source)
# eta <- matrix(c(0, 0, 0.2, 0.2, 0, 0.2, 0.2, 0, -0.2, 0.2), ncol = 2, byrow = TRUE)/10
gp <- c("homo xy", "heter xy")

X_pool_filtered <- read.csv('X_pool_filtered.csv')
B = 100# batch number
ITR = 2 #total iterations
t=1
eta <- matrix(c(0, 0, 0.1, 0.1, 0, 0.1, 0.1, 0, -0.1, 0.1), ncol = 2, byrow = TRUE)
n_train = 300
for (nonzero_beta2 in c(-2,-1,0,1,2)){
  for(t in c(3)){
    sim_result_10 = data.frame(matrix(ncol = 12, nrow = ITR*B))
    colnames(sim_result_10) <- c('p_s_sep_1','p_s_sep_2','p_s_sep','p_m_sep_1','p_m_sep_2','p_m_sep',
                                 'p_s_join_1','p_s_join_2','p_s_join','p_m_join_1','p_m_join_2','p_m_join')
    
    for(i in 1:ITR){
      Data_batch = DataGen_newhap(mc=B, n.train=300, n.test=3000, p=50, b0=1, 
                                  nonzero_beta1=1, nonzero_beta2=nonzero_beta2, gammas=eta[t, ], 
                                  var1=1, var2=1,var_D=1, seed=(i+1),group = "heter xy", X_pool=X_pool_filtered,snp_region = 1:50,eta_0 = 0.3)
      
      filename =  paste0('Data/Data_batch_0408_heter_datagen1_b2_',nonzero_beta2,'_', t,'setting', i,'.Rdata')
      #filename =  paste0('Data/Data_batch_0309_heter_beta1_eta_01_datagen1_nonzero_beta2_',nonzero_beta2,'_', t,'setting', i,'.Rdata')
      # filename =  paste0('Data/Data_batch_0301_homo_beta1_eta_01_datagen1_b0_',b,'_', t,'setting', i,'.Rdata')
      #filename =  paste0('Data/result_batch_0220_homo_beta1_eta_01_datagen1_n_train_',n_train,'_', t,'setting', i,'.csv')
      #filename = paste0('Data/Data_batch_0218_heter_eta_01_datagen1_nonzero_beta2_',nonzero_beta2,'_', t,'setting', i,'.RData')
      saveRDS(Data_batch,filename)
      #Data_batch = readRDS(filename)
      print('batch')
      print(i)
      sim_result_10 <- run_simulation(Data_batch, sim_result_10, i)
      
    }
    filename =   paste0('Data/result_batch_0408_heter_datagen1_b2_',nonzero_beta2,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0309_heter_beta1_eta_01_datagen1_nonzero_beta2_',nonzero_beta2,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0218_heter_eta_01_datagen1_nonzero_beta2_',nonzero_beta2,'_', t,'setting', i,'.csv')
    write.csv(sim_result_10, filename)
  }
}
