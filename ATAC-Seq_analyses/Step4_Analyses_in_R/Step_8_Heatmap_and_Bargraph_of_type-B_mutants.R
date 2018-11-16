#Load programs
library(dplyr)
library(ggplot2)
library(reshape2)

#Creating a matrix and heatmap of reads within increased DARs
Count = read.delim("ROOTS_UP_COUNTS_FILE.txt", header = TRUE, check.names = FALSE, as.is = TRUE)

Count$geneId = Roots_FC1.2_up$geneId
Count$annotation = Roots_FC1.2_up$annotation

Matrix = as.matrix(Count[, c(4:7)])
row.names(Matrix) = Count$geneId

#Row normalized heatmap
my_pallette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

heatmap_scale = heatmap.2(Matrix, 
          margins = c(12,9), 
          trace = "none", 
          scale = "row",
          col = my_pallette, 
          labRow = "")
          
#Bargraph of counts in annotations by genotype
Count_melted = melt(Count, id.vars = c("chr", "start", "end", "geneId", "annotation"))

Count_melted$treatment = c(rep("control", 666), rep("BA", 666), rep("control", 666), rep("BA", 666))

Count_melted_summ = ddply(Count_melted, c("variable"),
                          summarise, 
                          N = length(value),
                          mean = mean(value),
                          sd = sd(value),
                          se = sd / sqrt(N)
                          )

Count_melted_summ$treatment = c("control", "BA", "control", "BA")

Count_melted_summ$treatment = factor(Count_melted_summ$treatment, levels = 
                                       c("control", "BA"))

Count_melted_summ$name = c(rep("wild-type", 2), rep("type-B", 2))

Count_melted_summ$name = factor(Count_melted_summ$name, 
                                levels = c("wild-type", "type-B"))

#ANOVA on read counts
model = lm(value ~ variable, Count_melted)
anova(model)

#Turky post hoc test anova
TUKEY = HSD.test(model, "variable")
TUKEY

#Plot results
ggplot(Count_melted_summ, aes(x = treatment, y = mean, fill = name, ymax = 500)) + 
  geom_bar(position = "dodge", 
           stat="identity", 
           color = "black") + 
  scale_fill_manual(values = c("white", "green3")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,
                position=position_dodge(0.9)) + 
                  theme_classic() + 
  geom_text(aes(label = c("b", "a", "b", "b")), hjust=0.5, 
                vjust=-2, position=position_dodge(0.9)) + 
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) +
 theme(axis.ticks.x = element_blank(), 
       legend.position = "none", 
       axis.title = element_blank(),
       axis.text.x = element_blank())
