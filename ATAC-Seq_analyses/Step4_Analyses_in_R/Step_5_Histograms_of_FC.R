#Plotting histogram distribution of fold changes in roots
Roots_FC1.2$absolute_log2FC = abs(Roots_FC1.2$log2FC)
Shoots_FC1.2$absolute_log2FC = abs(Shoots_FC1.2$log2FC)

Roots_FC1.2[mapply(is.infinite, 
                       Roots_FC1.2)] <- NA
Shoots_FC1.2[mapply(is.infinite, 
                       Shoots_FC1.2)] <- NA

mean(Roots_FC1.2$absolute_log2FC, rm.na = TRUE)
mean(Shoots_FC1.2$absolute_log2FC, rm.na = TRUE)

ggplot(Roots_FC1.2, aes(x = absolute_log2FC)) +
  geom_histogram(binwidth = 0.05, 
                 colour= "black", 
                 fill = "red") + 
  theme_classic() + 
  geom_vline(aes(xintercept = mean(absolute_log2FC, 
                                   rm.na = TRUE)), 
             color = "black", linetype = "dashed", 
             size = 1)

ggplot(Shoots_FC1.2, aes(x = absolute_log2FC)) +
  geom_histogram(binwidth = 0.05, 
                 colour= "black", 
                 fill = "red") + 
  theme_classic() + 
  geom_vline(aes(xintercept = mean(absolute_log2FC, 
                                   rm.na = TRUE)), 
             color = "black", linetype = "dashed", 
             size = 1)  
