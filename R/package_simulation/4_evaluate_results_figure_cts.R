library(ggplot2)
i=100
t=3
i=1
t=5
t=1
i=34
i=2
#n_train in c(150, 200, 250,300)
for (nonzero_beta2 in c(-2,-1,0,1,2)){
  for(t in c(3)){
    #print(b)
    print(t) 
    filename =   paste0('Data/result_batch_0408_heter_datagen1_b2_',nonzero_beta2,'_', t,'setting', i,'.csv')
    #filename = paste0('Data/result_batch_0313_heter_binom_nonzero_beta2_',nonzero_beta2,'_', t,'setting',i,'.csv')
    #filename =  paste0('Data/result_batch_0309_heter_beta1_eta_01_datagen1_nonzero_beta2_',nonzero_beta2,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0307_heter_beta1_eta_01_datagen1_nonzero_beta1_',nonzero_beta1,'_', t,'setting', i,'.csv')
    #filename = paste0('Data/result_batch_0301_homo_beta1_eta_01_datagen1_b0_',b,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0307_b12b01_heter_test_',b,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0218_heter_beta1_eta_01_datagen1_nonzero_beta1_',nonzero_beta1,'_', t,'setting', i,'.csv')
    # filename =  paste0('Data/result_batch_0216_heter_beta1_eta_01_datagen1_b0_',b,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0220_homo_beta1_eta_01_datagen1_n_train_',n_train,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0216_homo_beta1_eta_01_datagen1_b0_',b,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0212_heter_beta12_datagen1_b0_',b,'_', t,'setting', i,'.csv')
    # filename =  paste0('Data/result_batch_0210_eta0_03_homo_beta12_datagen1_ntrain_',n_train,'_', t,'setting', i,'.csv')
    # filename =  paste0('Data/result_batch_0210_eta0_03_beta12_datagen1_ntrain_',n_train,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta12_datagen1_ntrain_',n_train,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_homo_beta12_datagen1_ntrain_',n_train,'_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta12_datagen1_rb23',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta11_datagen1_rb22',t,'setting', i,'.csv')
    # filename =  paste0('Data/result_batch_0210_eta02_heter_beta13_datagen1_rb22',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta13_datagen1_rb22',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta13_datagen1',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0210_eta02_heter_beta02_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0205_eta02_heter_beta22_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0205_eta005_heter_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0205_eta005_heter_pi23_',n_train,'_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0205_eta005_heter_var001_',n_train,'_',t,'setting', i,'.csv')
    #filename = paste0('Data/result_batch_0205_eta005_heter_b_',b, '_', t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_b01_eta02_homo_3000_b_',b,'_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_b01_eta02_heter_3000_b_',b,'_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_b01_eta02_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_b01_eta02_homo_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_beta1_eta02_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_beta02_b02_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_varD2_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_eta5_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_varD10_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_eta03_1000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_beta01_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/result_batch_0203_newhap_center_varD10_3000_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0202_newhap_nonenter_eta03_beta2_2_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0128_newhap_center_eta03_beta2_1_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0202_newhap_center_eta03_beta2_3_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0202_newhap_center_eta03_beta2_2_',t,'setting', i,'.csv')
     #filename =  paste0('Data/Data_batch_0121_newhap_center_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0128_newhap_center_',t,'setting', i,'.csv')
    # filename =  paste0('Data/Data_batch_0128_newhap_center_eta03_var01_',t,'setting', i,'.csv')
    #filename =  paste0('Data/Data_batch_0128_newhap_center_eta03_',t,'setting', i,'.csv')
    print(filename)
   
    s_2<-read.csv(filename, row.names = NULL)[, -1]
    s_3 <- s_2[which(!is.na(s_2$p_m_join) & !is.na(s_2$p_s_join) & !is.na(s_2$p_s_join) & !is.na(s_2$p_s_sep)), ]
    #s_3 <-s_3[1:4000, ]
    # print(head(s_3))
    s_4 <- s_3[which(!is.infinite(-log10(s_3$p_s_join)) & !is.infinite(-log10(s_3$p_m_join))), ]
    s_join  = length(which(s_3$p_s_join < 0.05))/nrow(s_3)
   # print(s_join)
    m_join  = length(which(s_3$p_m_join < 0.05))/nrow(s_3)
    #print(m_join)
    s_join  = length(which(s_3$p_s_join_1 < 0.05))/nrow(s_3)
    #print(s_join)
    m_join  = length(which(s_3$p_s_join_2 < 0.05))/nrow(s_3)
    #print(m_join)
    s_sep  = length(which(s_3$p_s_sep < 0.05))/nrow(s_3)
   # print(s_sep)
    m_sep  = length(which(s_3$p_m_sep < 0.05))/nrow(s_3)
    #print(m_sep)
    
    # s_3 <- s_2[which(!(is.infinite(-log10(s_2$p_m_join)) & is.infinite(-log10(s_2$p_s_join)))), ]
    #c <- cor(s_3$p_m_join, s_3$p_s_join, use="complete.obs")
    #print(c)
    #c <- cor(-log(s_3$p_m_join+1e-200), -log(s_3$p_s_join+1e-200), use="complete.obs")
    #print(c)
    #c <- cor(-log10(s_3$p_m_join), -log10(s_3$p_s_join), use="complete.obs")
    c <- cor(-log10(s_4$p_m_sep), -log10(s_4$p_m_join), use="complete.obs")
    print(c)
    c <- cor(-log10(s_4$p_s_sep), -log10(s_4$p_s_join), use="complete.obs")
    print(c)
    c <- cor(-log10(s_4$p_m_sep), -log10(s_4$p_s_sep), use="complete.obs")
    print(c)
    c <- cor(-log10(s_4$p_m_join), -log10(s_4$p_s_join), use="complete.obs")
    print(c)
  }
}
    # Prepare data frame with -log10 transformation
    df_qq <- data.frame(
      Mixcan = -log10(sort(s_4$p_m_join)),
      SMixcan = -log10(sort(s_4$p_s_join))
    )
    
    # Compute correlation coefficient
    cor_value <- cor(df_qq$Mixcan, df_qq$SMixcan, method = "pearson")
    cor_text <- paste0("Pearson r = ", round(c, 3))
    
    # Create QQ plot
    ggplot(df_qq, aes(x = Mixcan, y = SMixcan)) +
      geom_point(alpha = 0.6, color = "blue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Reference line y = x
     # Set y-axis limit
      theme_minimal() +
      labs(
           x = "MiXcan -log10(p)",
           y = "S-MiXcan -log10(p)") +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)) +
      annotate("text", x = 15, y = 0.5, label = cor_text, size = 5, hjust = 1, color = "black")  # Display correlation
  }
}


# Function to prepare QQ plot data
qq_data <- function(p_values, method_name) {
  data.frame(
    observed = sort(-log10(p_values)),  # Observed -log10(p)
    expected = sort(-log10(ppoints(length(p_values)))),  # Expected -log10(p)
    method = method_name
  )
}

# Prepare data for both sets
df_qq1 <- qq_data(p_values_1, "Method 1")
df_qq2 <- qq_data(p_values_2, "Method 2")

# Combine into a single dataframe
df_qq <- rbind(df_qq1, df_qq2)

# Create QQ plot
ggplot(df_qq, aes(x = expected, y = observed, color = method)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "QQ Plot of P-values",
       x = "Expected -log10(p)",
       y = "Observed -log10(p)") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

b <- c(-2,-1,0,1,2)
B=100
ITR=1
n_train=300
n_test=3000
beta1=1
beta2=1
grp = 'heter'




  result_summary <- data.frame(matrix(ncol = 11, nrow = 5*11))
  colnames <- c('n_replicates', 'n_train', 'n_test','b0','nonzero_beta1','nonzero_beta2','group','gamma','power_S','power_M','cor')
  for(t in 1:5){
  filename =  paste0('Data/result_batch_0203_newhap_center_b01_eta02_homo_3000_b_',b,'_',t,'setting', i,'.csv')
  result_summary[, 'n_replicates'] <- ITR*B
  result_summary[, 'n_train'] <- n_train
  result_summary[, 'n_test'] <- n_test
  result_summary[t, 'b0'] <- 1
  result_summary[, 'nonzero_beta1'] <- beta1
  result_summary[, 'nonzero_beta2'] <- beta2
  result_summary[, 'group'] <- grp
  result_summary[, 'gamma'] <- str(eta[t, ])
  
  s_2<-read.csv(filename, row.names = NULL)[, -1]
  s_3 <- s_2[which(!is.na(s_2$p_m_join) & !is.na(s_2$p_s_join) & !is.na(s_2$p_s_join) & !is.na(s_2$p_s_sep)), ]
  s_4 <- s_3[which(!is.infinite(-log10(s_3$p_s_join)) & !is.infinite(-log10(s_3$p_m_join))), ]
  s_join  = length(which(s_3$p_s_join < 0.05))/nrow(s_3)
  print(s_join)
  m_join  = length(which(s_3$p_m_join < 0.05))/nrow(s_3)
  print(m_join)
  result_summary[t, 'power_S'] <- s_join
  result_summary[t, 'power_M'] <- m_join
  result_summary[t, 'cor'] <- cor(-log10(s_4$p_m_join), -log10(s_4$p_s_join), use="complete.obs")
  }
  
  filename_out =  paste0('Data/summary_0203_heter_var001','.csv')
  write.csv(result_summary, filename_out)



i=2
t=1
t=5
for(t in 1:5){
  filename =  paste0('Data/result_batch_0203_newhap_center_varD10_3000_',t,'setting', i,'.csv')
  s_2<-read.csv(filename)
  s_3 <- s_2[which(!is.na(s_2$p_m_join)), ]
  s_4 <- s_3[which(!is.infinite(-log10(s_3$p_s_join)) & !is.infinite(-log10(s_3$p_m_join))), ]
  # s_3 <- s_2[which(!(is.infinite(-log10(s_2$p_m_join)) & is.infinite(-log10(s_2$p_s_join)))), ]
  c <- cor(-log10(s_4$p_m_join), -log10(s_4$p_s_join), use="complete.obs")
  print(c)
}


plot(-log10(sim_result_10$p_m_join), -log10(sim_result_10$p_s_join), main='cts')
abline(0,1)

sim_result_pw <- sim_result_10
cor(-log10(sim_result_10$p_m_join), -log10(sim_result_10$p_s_join), use="complete.obs")

length(which(sim_result_10$p_s_join < 0.05))
length(which(sim_result_10$p_m_join < 0.05))

sim_result_bin_s <- sim_result_10[which(!is.na(sim_result_10$p_s_sep_1)), ]

s_join = length(which(sim_result_bin_s$p_s_join < 0.05))/nrow(sim_result_10)
s_join_1  = length(which(sim_result_bin_s$p_m_join < 0.05))/nrow(sim_result_10)
length(which(sim_result_bin_s$p_s_join_1 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_m_join_1 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_s_join_2 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_m_join_2 < 0.05))/nrow(sim_result_bin_s)

s_join  = length(which(sim_result_bin_s$p_s_join < 0.05))/nrow(sim_result_10)
s_join_1 = length(which(sim_result_bin_s$p_s_join_1 < 0.05))/nrow(sim_result_10)
s_join_2 = length(which(sim_result_bin_s$p_s_join_2 < 0.05))/nrow(sim_result_10)

s_sep = length(which(sim_result_bin_s$p_s_sep < 0.05))/nrow(sim_result_10)
s_sep_1 = length(which(sim_result_bin_s$p_s_sep_1 < 0.05))/nrow(sim_result_10)
s_sep_2 = length(which(sim_result_bin_s$p_s_sep_2 < 0.05))/nrow(sim_result_10)

barplot(height = c(s_join, s_join_1, s_join_2, s_sep, s_sep_1, s_sep_2), names.arg = c('s_join', 's_join_1','s_join_2', 's_sep', 's_sep_1', 's_sep_2'), cex.names = 0.8, ylim=c(0,0.05), main='power')
