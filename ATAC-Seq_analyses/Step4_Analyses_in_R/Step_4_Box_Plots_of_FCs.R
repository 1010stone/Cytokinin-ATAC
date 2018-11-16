#Plotting box plots of fold changes in roots
#Change exons to display just exon
Roots_FC1.2$annotation[grep(pattern = "Exon", ignore.case = T, x = Roots_FC1.2$annotation)] = "Other Exon"
Shoots_FC1.2$annotation[grep(pattern = "Exon", ignore.case = T, x = Shoots_FC1.2$annotation)] = "Other Exon"

Roots_FC1.2$annotation[grep(pattern = "Downstream", ignore.case = T, x = Roots_FC1.2$annotation)] = "Downstream (<=3kb)"
Shoots_FC1.2$annotation[grep(pattern = "Downstream", ignore.case = T, x = Shoots_FC1.2$annotation)] = "Downstream (<=3kb)"

Roots_FC1.2$annotation = 
  factor(Roots_FC1.2$annotation, 
         c("Promoter (<=1kb)", 
          "Promoter (1-2kb)", 
          "Promoter (2-3kb)", 
          "3' UTR", 
          "Other Exon",
          "Downstream (<=3kb)", 
          "Distal Intergenic"))
Shoots_FC1.2$annotation = 
  factor(Shoots_FC1.2$annotation, 
         c("Promoter (<=1kb)", 
          "Promoter (1-2kb)", 
          "Promoter (2-3kb)", 
          "3' UTR", 
          "Other Exon",
          "Downstream (<=3kb)", 
          "Distal Intergenic"))

ggplot(Roots_FC1.2, aes(x = annotation,
                            y = absolute_log2FC, 
                            fill = annotation)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  scale_fill_manual(values=c("#a6cee3", 
                             "#b2df8a", 
                             "#CDAF95", 
                             "#EE4000",
                             "#ff7f00",
                             "#6a3d9a",
                             "#b15928"))



ggplot(Shoots_FC1.2, aes(x = annotation,
                            y = absolute_log2FC, 
                            fill = annotation)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  scale_fill_manual(values=c("#a6cee3", 
                             "#b2df8a", 
                             "#CDAF95", 
                             "#EE4000",
                             "#ff7f00",
                             "#6a3d9a",
                             "#b15928"))

