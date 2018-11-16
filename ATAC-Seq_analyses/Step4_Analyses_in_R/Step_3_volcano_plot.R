#Volcano plot of showing all differential regoins from merged datasets

#Highlight genes with FC greater or equal to 1.2 and adjusted p-value less than or equal to 0.05

Roots$threshold = as.factor(abs(Roots$log2FC) >= 0.2630344 & abs(Roots$padj) <= 0.05)
Shoots$threshold = as.factor(abs(Shoots$log2FC) >= 0.2630344 & abs(Shoots$padj) <= 0.05)

#Construct the volcano plots

ggplot(data = Roots, aes(x = log2FC, y = -log10(pval), colour = threshold)) +
  geom_point(alpha = 0.4, size = 1) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() + theme(legend.position = "none")

ggplot(data = Shoots, aes(x = log2FC, y = -log10(pval), colour = threshold)) +
  geom_point(alpha = 0.4, size = 1) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() + theme(legend.position = "none")
