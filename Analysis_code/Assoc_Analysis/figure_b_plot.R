library(bacon)

# Directory where your files are stored
path <- "/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/DRIVE/"


# Create list of file paths
file_list <- paste0(path, "drive_result_chr", sprintf("%02d", 3:11), "_recover.csv")

# Read and combine all files into one data frame
all_results <- do.call(rbind, lapply(file_list, read.csv))

names(all_results)

hist(all_results$p_m_join, breaks = 50, main = "Joint P-value of MiXcan in DRIVE", xlab = "P-value")

pvals <- all_results$p_s_join
expected <- -log10(ppoints(length(pvals)))
observed <- -log10(sort(pvals))

# figure c--------------------
plot(expected, observed,
     main = "QQ Plot of joint SMiXcan P-values in DRIVE",
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
     pch = 20, cex = 0.6)
abline(0, 1, col = "red", lty = 2)

p_sorted <- sort(pvals)
n <- length(p_sorted)
p <- all_results$p_s_join
y=qnorm(1 - p)
bc <- bacon(y)
estimates(bc)
inflation(bc)
newp=pval(bc)

# change page size






