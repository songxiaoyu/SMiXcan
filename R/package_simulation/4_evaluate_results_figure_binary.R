

xt = sim_data$x.train
MAF <- colSums(X_pool)/(nrow(X_pool)*2)
X_pool_filtered <- X_pool[, MAF > 0.01]
View(X_pool)




pw_result <- sim_result_bin


plot(-log10(sim_result_bin$p_m_join), -log10(sim_result_bin$p_s_join), xlab = '-log10(p_mixcan)', ylab='-log10(p_s-mixcan)',main='Separate Analysis (tissue level)')
abline(0,1)
cor(-log10(sim_result_bin$p_m_join), -log10(sim_result_bin$p_s_join), use="complete.obs")
text(3,35,'cor = 0.999')
text(1,8,'cor = 0.997')

s_join  = length(which(pw_result$p_s_join < 0.05))/nrow(pw_result)
s_join_1 = length(which(pw_result$p_s_join_1 < 0.05))/nrow(pw_result)
s_join_2 = length(which(pw_result$p_s_join_2 < 0.05))/nrow(pw_result)

s_sep = length(which(pw_result$p_s_sep < 0.05))/nrow(pw_result)
s_sep_1 = length(which(pw_result$p_s_sep_1 < 0.05))/nrow(pw_result)
s_sep_2 = length(which(pw_result$p_s_sep_2 < 0.05))/nrow(pw_result)

barplot(height = c(s_join, s_join_1, s_join_2, s_sep, s_sep_1, s_sep_2), names.arg = c('s_join', 's_join_1','s_join_2', 's_sep', 's_sep_1', 's_sep_2'), cex.names = 0.8, ylim=c(0,1), main='power (S-MiXcan)')

m_join  = length(which(pw_result$p_m_join < 0.05))/nrow(pw_result)
m_join_1 = length(which(pw_result$p_m_join_1 < 0.05))/nrow(pw_result)
m_join_2 = length(which(pw_result$p_m_join_2 < 0.05))/nrow(pw_result)

s_sep = length(which(pw_result$p_s_sep < 0.05))/nrow(pw_result)
s_sep_1 = length(which(pw_result$p_s_sep_1 < 0.05))/nrow(pw_result)
s_sep_2 = length(which(pw_result$p_s_sep_2 < 0.05))/nrow(pw_result)

barplot(height = c(s_join, s_join_1, s_join_2, m_join, m_join_1, m_join_2), names.arg = c('s_join', 's_join_1','s_join_2', 'm_join', 'm_join_1','m_join_2'), cex.names = 0.8, ylim=c(0,1), main='power')

pw_result
s_1 <- read.csv(paste0('Data/result_0121_newhap_t1e_binom_center',2,'.csv'))
sim_result_bin_s <- sim_result_bin[which(!is.na(s_1$p_s_sep_1)), ]
sim_result_bin_s <- sim_result_bin_s[1:2000, ]

s_join = length(which(sim_result_bin_s$p_s_join < 0.05))/nrow(sim_result_bin_s)
s_join_1  = length(which(sim_result_bin_s$p_m_join < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_s_join_1 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_m_join_1 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_s_join_2 < 0.05))/nrow(sim_result_bin_s)
length(which(sim_result_bin_s$p_m_join_2 < 0.05))/nrow(sim_result_bin_s)

s_join  = length(which(sim_result_bin_s$p_s_join < 0.05))/nrow(sim_result_bin_s)
s_join_1 = length(which(sim_result_bin_s$p_s_join_1 < 0.05))/nrow(sim_result_bin_s)
s_join_2 = length(which(sim_result_bin_s$p_s_join_2 < 0.05))/nrow(sim_result_bin_s)

s_sep = length(which(sim_result_bin_s$p_s_sep < 0.05))/nrow(sim_result_bin_s)
s_sep_1 = length(which(sim_result_bin_s$p_s_sep_1 < 0.05))/nrow(sim_result_bin_s)
s_sep_2 = length(which(sim_result_bin_s$p_s_sep_2 < 0.05))/nrow(sim_result_bin_s)

barplot(height = c(s_join, s_join_1, s_join_2, m_join, m_join_1, m_join_2), names.arg = c('s_join', 's_join_1','s_join_2', 's_sep', 's_sep_1', 's_sep_2'), cex.names = 0.8, ylim=c(0,0.05), main='power')

m_join  = length(which(sim_result_bin_s$p_m_join < 0.05))/nrow(sim_result_bin_s)
m_join_1 = length(which(sim_result_bin_s$p_m_join_1 < 0.05))/nrow(sim_result_bin_s)
m_join_2 = length(which(sim_result_bin_s$p_m_join_2 < 0.05))/nrow(sim_result_bin_s)

m_sep = length(which(s_1$p_m_sep < 0.05))/nrow(s_1)
m_sep_1 = length(which(s_1$p_m_sep_1 < 0.05))/nrow(s_1)
m_sep_2 = length(which(s_1$p_m_sep_2 < 0.05))/nrow(s_1)
barplot(height = c(m_join, m_join_1, m_join_2, m_join, m_join_1, m_join_2), names.arg = c('m_join', 'm_join_1','m_join_2', 'm_sep', 'm_sep_1', 'm_sep_2'), cex.names = 0.8, ylim=c(0,0.05), main='power')


barplot(height = c(m_join, m_sep, s_join, s_sep), names.arg = c('m_join', 'm_sep', 's_join', 's_sep'), cex.names = 0.8, ylim=c(0,1), main='Type I error')
bar_heights <- c(m_join, m_sep, s_join, s_sep)  # Heights of the bars
bar_positions <- barplot(
  height = bar_heights,
  names.arg = c('m_join', 'm_sep', 's_join', 's_sep'),
  cex.names = 0.8,
  ylim = c(0, 1),
  main = 'Type I Error',
  ylab = 'Type I Error',
)

# Add values on top of the bars
text(
  x = bar_positions,            # x-coordinates of bars
  y = bar_heights,              # y-coordinates of bars
  labels = round(bar_heights-0.0001, 3), # Labels (rounded to 2 decimal places)
  pos = 3,                      # Position above the bars
  cex = 0.6,                    # Text size
  col = "black"                 # Text color
)
abline(h = 0.05, col = "red", lty = 2, lwd = 2)  # `h` specifies the y-value


setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_Simulation')

s_1 <- read.csv(paste0('Data/result_1230_newhap_t1e_binom_center',20,'.csv'))
# i=5

i=2
for(t in 1:5){
  s_2<-read.csv(paste0('Data/result_0121_newhap_binom_',t,'setting',i,'.csv'))
  s_3 <- s_2[which(!is.na(s_2$p_m_join)), ]
  s_4 <- s_3[which(!is.infinite(-log10(s_3$p_s_join))), ]
  # s_3 <- s_2[which(!(is.infinite(-log10(s_2$p_m_join)) & is.infinite(-log10(s_2$p_s_join)))), ]
  c <- cor(-log10(s_4$p_m_join), -log10(s_4$p_s_join), use="complete.obs")
  print(c)
}

plot(-log10(s_4$p_m_join), -log10(s_4$p_s_join),xlim=c(0,200),ylim=c(0,200))

