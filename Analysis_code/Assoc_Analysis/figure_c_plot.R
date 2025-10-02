
bcac_result<- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/BCAC/bcac2020_result/combined_with_cytoband.csv")
head(bcac_results)
# Extract p-values
pvals <- bcac_result$p_s_join

# Expected vs Observed
expected <- -log10(ppoints(length(pvals)))
observed <- -log10(sort(pvals))

# run genome inflation 
# Genomic inflation factor
lambda <- 1.007253

# Save
#pdf("QQplot of .pdf", width=5, height=5)

# QQ plot
plot(expected, observed,
     xlab=expression(Expected~~-log[10](italic(p))),
     ylab=expression(Observed~~-log[10](italic(p))),
     main="QQ Plot of SMiXcan",
     pch=19, cex=0.6, col="blue", las=1)
abline(0, 1, col="red", lwd=2, lty=2)

# Add lambda to the plot (top-left corner)
text(x=min(expected), y=max(observed)*0.9,
     labels=bquote(lambda["GC"] == .(round(lambda,3))),
     adj=0, cex=1)

dev.off()

hist(bcac_result$p_s_join, breaks=50,
     main="Histogram of p-values", xlab="p-value", col="skyblue")



