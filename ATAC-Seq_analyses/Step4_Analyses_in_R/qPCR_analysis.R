# ANOVA and Tukey post hoc analysis
library(agricolae)
library(userfriendlyscience)
library(car)
library(DTK)
library(DescTools)
library(multcomp)
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)

#Load in qPCR data
Data = read.csv("ATAC-qPCR/20170606 - ATAC-qPCR summary.csv", header = TRUE)

#Calculate means and perform ANOVA
index = factor(Data$Samples, levels = c("R-", "R+", "B-", "B+"))
tapply(X = Data$FCs, INDEX = index, FUN = mean)
model = lm(FCs ~ Samples, Data)
anova(model)
summary(model)

#muliple comparison using multcomp package
mc = glht(model, mcp(Samples = "Tukey"))

#for Tukey posthoc analysis, then just summarize mc
summary(mc)

cld(mc, level = 0.05, decreasing = TRUE)

Data_Summary = Data %>%
  group_by(Samples, Treatment) %>%
  summarise(mean_FCs = mean(FCs), 
            sd_FCs = sd(FCs), 
            n_FCs = n(),
            se_FCs = sd(FCs)/sqrt(n()))

Data_Summary$Samples = factor(Data_Summary$Samples, 
                                 levels = c("R-", "R+", "B-", "B+"))

Data_Summary$Treatment = factor(Data_Summary$Treatment, 
                                 levels = c("control", "BA"))

Data_Summary$letter = factor(c("c", "bc", "b", "a"), levels = c("c", "bc", "b", "a"))

Data$Treatment = factor(Data$Treatment, 
                                 levels = c("control", "BA"))

Data$Samples = factor(Data$Samples, 
                         levels = c("R-", "R+", "B-", "B+"))

ggplot(Data_Summary, 
       aes(x = Treatment, 
           y = mean_FCs, 
           fill = Samples,
           ymax = 2.25, ymin = 0)) +
  geom_bar(position = "dodge", 
           stat="identity", 
           color = "black") +
  scale_fill_manual(values = c("white", "white", "green3", "green3")) +
  geom_text(aes(label = letter), 
           hjust=0.5,
           vjust=-4, position=position_dodge(0.9)) +
  geom_point(data = Data, 
             aes(x = Treatment, y = FCs, fill = Samples), 
             position = position_dodge(0.9), 
             stat = "identity", 
         #    pch = 21, 
              size = 2, 
         #    stroke = 1, 
         #    alpha = 0.7
             color = "black") +
  geom_point(data = Data, 
             aes(x = Treatment, y = FCs, fill = Samples), 
             position = position_dodge(0.9), 
             stat = "identity", 
         #    pch = 21, 
              size = 1, 
         #    stroke = 1, 
         #    alpha = 0.7
             color = "grey") +
  geom_errorbar(aes(ymax=mean_FCs+se_FCs, 
                    ymin=mean_FCs-se_FCs),
                width = 0.2, size = 0.5,
                color = "black",
                position=position_dodge(0.9)) +
 theme_classic() +
 scale_y_continuous(expand = c(0,0)) +
 theme(axis.ticks.x = element_blank(), 
       legend.position = "none", 
       axis.title = element_blank(),
       axis.text = element_blank())
