# This script compute and compare overall methylation statistics

library(methylKit)
library(genomation)
library(stringr)
library(Biobase)
library(pheatmap)
library(tidyr)
source("read_bismark_methyl.R")
source("meth_analysis_funcs.R")

strains <- factor(c(rep("EH1516",2), rep("EH217", 3)), levels=c("EH1516", "EH217"))

# Calculate p-value for t-test the mean methylation levels between EH1516 and EH217
calc_avg_meth_pval <- function(meth.obj, n_permutations = 10000) {
  cov.cols <- c("coverage1", "coverage2", "coverage3", "coverage4", "coverage5")
  numTs.cols <- c("numCs1", "numCs2", "numCs3", "numCs4", "numCs5")
  # Calculate the methylation levels
  meth.level <- getData(meth.obj)[, numTs.cols] / getData(meth.obj)[, cov.cols]
  #cat(colMeans(meth.level))
  meth.level.EH1516 <- meth.level[, strains == "EH1516"]
  meth.level.EH217 <- meth.level[, strains == "EH217"]
  
  group1 <- as.vector(rowMeans(meth.level.EH1516))
  group2 <- as.vector(rowMeans(meth.level.EH217))
  
  # Should not use t-test because methylation isn't in normal distribution
  #res.t.test <- t.test(group1, group2)
  res.wilcox <- wilcox.test(group1, group2, paired=TRUE)
  
  # Calculate the observed difference in means
  # obs_diff <- mean(group1) - mean(group2)
  # 
  # # Combine the data
  # combined_data <- c(group1, group2)
  # perm_diffs <- numeric(n_permutations)
  # 
  # # Perform the permutations
  # for (i in 1:n_permutations) {
  #   print(paste("permutation: ", i))
  #   permuted_data <- sample(combined_data)  # Shuffle the combined data
  #   perm_group1 <- permuted_data[1:length(group1)]  # First half is new group 1
  #   perm_group2 <- permuted_data[(length(group1)+1):length(combined_data)]  # Second half is new group 2
  #   
  #   # Calculate the difference in means for the permuted groups
  #   perm_diffs[i] <- mean(perm_group1) - mean(perm_group2)
  # }
  # 
  # # Calculate p-value: proportion of permuted differences greater than or equal to observed difference
  # p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  res <- list(means = colMeans(meth.level), pval=res.wilcox)
  return (res)
}
# meth.CHG.cov10 has already been made available previously
# Test statistical significance of EH1516 and EH217 have different mean
# permutation test results in pval = 0
res.CHG <- calc_avg_meth_pval(meth.CHG.cov10)
print(res.CHG$means)
print(res.CHG$pval)

res.CHG.wilcox <- calc_avg_meth_pval(meth.CHG.cov10)
print(res.CHG.wilcox$means)
print(res.CHG.wilcox$pval)

# test for CpG methylation using permutation test
# resulting pval = 0
res.CpG <- calc_avg_meth_pval(meth10)
print(res.CpG$means)
print(res.CpG$pval)

res.CpG.wilcox <- calc_avg_meth_pval(meth10)
print(res.CpG.wilcox$means)
print(res.CpG.wilcox$pval)

# Plot bar chart of methylation levels in various contexts
# meth.Mean <- cbind(res.CpG.wilcox$means, res.CHG.wilcox$means)
# row.names(meth.Mean) <- c("EH1516B", "EH1516C", "EH217A", "EH217B", "EH217C")
# colnames(meth.Mean) <- c("CpG", "CHG")
# meth.Mean <- as.data.frame(meth.Mean)
# meth.Mean$strain <- strains #c("EH1516", "EH1516", "EH217", "EH217", "EH217")

# library(reshape2)
# library(ggplot2)
# mdat <- melt(meth.Mean, id.vars = "strain")
# df_long <- df %>% 
#   gather(key = "Type", value = "Value", CpG, CHG)
# #rm(meth.CHG.EH1516.vec, meth.CHG.EH217)
# ggplot(mdat, aes(variable, value, fill=strain)) + 
#   geom_bar(stat="identity", position="dodge")
# color_mapping <- c("EH1516" = "red", "EH217"="blue")
# bar_color <- color_mapping[meth.Mean$strain]
# barplot(as.matrix(meth.Mean[, c(1, 2)]), beside=TRUE, col=bar_color)

