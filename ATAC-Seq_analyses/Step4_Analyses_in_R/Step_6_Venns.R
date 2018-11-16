# Load programs
library(xlsx)
library(GenomicFeatures)
library(ChIPseeker)
library(gdata)
library(gplots)
library(VennDiagram)
library(ggplot2)

#Load in datasets
tair10 <- loadDb("TAIR10_GFF3_genes_chr_removed.sqlite")
RootsGenes = read.delim("RootsGenes.txt", header = T, 
                        col.names = c("logFC", "logCPM", "LR",
                                      "PValue", "FDR", "gene"))

ShootsGenes = read.delim("ShootsGenes.txt", header = T, 
                         col.names = c("logFC", "logCPM", 
                                       "LR", "PValue", "FDR",
                                       "gene"))

#Separate induced/repressed genes from RNA-Seq
RootsGenes_up = RootsGenes[RootsGenes$logFC > 0.5, ]
RootsGenes_down = RootsGenes[RootsGenes$logFC < 0.5, ]

ShootsGenes_up = ShootsGenes[ShootsGenes$logFC > 0.5, ]
ShootsGenes_down = ShootsGenes[ShootsGenes$logFC < 0.5, ]

#Load in differential ATAC regions and genes from roots
Roots_all_FC1.2 = readPeakFile("Roots_WT_All_diffRepsFC1.2.bed")
Roots_all_FC1.2_up = readPeakFile("Roots_WT_UP_diffRepsFC1.2.bed")
Roots_all_FC1.2_down = readPeakFile("Roots_WT_DOWN_diffRepsFC1.2.bed")

Roots_ATAC_genes_all = read.xlsx("Roots_WT_All_diffrepsFC1.2.xlsx", header = T, 
                            sheetIndex = 1)
Roots_ATAC_genes_UP = Roots_ATAC_genes_all[Roots_ATAC_genes_all$FC > 1.2,]
Roots_ATAC_genes_DOWN = Roots_ATAC_genes_all[Roots_ATAC_genes_all$FC < 0.8, ]

#Load in differential ATAC regions and genes from shoots
Shoots_all_FC1.2 = readPeakFile("Shoots_WT_All_diffRepsFC1.2.bed")
Shoots_all_FC1.2_up = readPeakFile("Shoots_WT_UP_diffRepsFC1.2.bed")
Shoots_all_FC1.2_down = readPeakFile("Shoots_WT_DOWN_diffRepsFC1.2.bed")

Shoots_ATAC_genes_all = read.xlsx("Shoots_WT_All_diffrepsFC1.2.xlsx", header = T, 
                            sheetIndex = 1)
Shoots_ATAC_genes_UP = Shoots_ATAC_genes_all[Shoots_ATAC_genes_all$FC > 1.2,]
Shoots_ATAC_genes_DOWN = Shoots_ATAC_genes_all[Shoots_ATAC_genes_all$FC < 0.8, ]

#Load in type-B ARR ChIP files
arr1_GR = readPeakFile("GSM2476320_ARR1.2ypet_LD_6BA.peak.bed")
arr10_GR = readPeakFile("arr10_overlaps.bed") #generated using bedops intersect of Zubo and Ecker arr10 files
arr12_GR = readPeakFile("GSM2476326_ARR12Ypet_6BA.peak.bed")

#Venn of Roots ATAC and RNA genes
up = list("Roots_ATAC_up" = Roots_ATAC_genes_UP$geneId, 
          "Roots_RNA_up" = RootsGenes_up$gene)
down = list("Roots_ATAC_down" = Roots_ATAC_genes_DOWN$geneId, 
            "Roots_RNA_down" = RootsGenes_down$gene)
venn_up = venn.diagram(up, cex = 1, col = c("red", "black"), ext.percent = 0, 
                       filename = NULL)
venn_down = venn.diagram(down, cex = 1, col = c("red", "black"), ext.percent = 0, 
                         filename = NULL)
grid.draw(venn_up)
phyper(157, 471, 25000-471, 471+585, lower.tail = FALSE, log.p = F)
grid.draw(venn_down)
phyper(20, 20+59, 25000-20-59, 374+20+59, lower.tail = FALSE, log.p = F)

#Venn of Shoots ATAC and RNA genes
up = list("Shoots_RNA_up" = ShootsGenes_up$gene, "Shoots_ATAC_up" = Shoots_ATAC_genes_UP$geneId)
down = list("Shoots_ATAC_down" = Shoots_ATAC_genes_DOWN$geneId, 
            "Shoots_RNA_down" = ShootsGenes_down$gene)
venn_up = venn.diagram(up, cex = 1, col = c("red", "black"), ext.percent = 0, 
                       filename = NULL)
venn_down = venn.diagram(down, cex = 1, col = c("red", "black"), ext.percent = 0, 
                         filename = NULL)
grid.draw(venn_up)
grid.draw(venn_down)
phyper(7, 56+7, 25000-56-7, 139+7+56, lower.tail = FALSE, log.p = F)

#Overlap of root and shoot differentially expressed genes 
#Venn of overlaps
Gene_Overlaps = list("Roots" = RootsGenes$gene, 
                     "Shoots" = ShootsGenes$gene)
venn_overlap = venn.diagram(Gene_Overlaps, cex = 1, col = c("red", "black"), 
                            ext.percent = 0, filename = NULL)
grid.draw(venn_overlap)

#Overlap of root and shoot differential ATAC regions 
#Venn of overlap
a = sum(Roots_all_FC1.2 %over% Shoots_all_FC1.2 =="FALSE")
b = sum(Shoots_all_FC1.2 %over% Roots_all_FC1.2 =="FALSE")
overlaps = numOverlaps(Roots_all_FC1.2, Shoots_all_FC1.2)

venn <- draw.pairwise.venn(area1 = a + overlaps, area2 = b + overlaps, 
                           cross.area = overlaps, cex = 1, 
                           col = c("darkmagenta", "darkgoldenrod"), ext.percent = 0)
grid.draw(venn)

#Venn of type-B ARR datasets with DARs
overlaps_roots = list("arr1" = arr1_GR, "arr10" = arr10_GR, 
                      "arr12" = arr12_GR, "Roots" = Roots_all_FC1.2)
vennplot(overlaps_roots, by = "gplots")

overlaps_shoots = list("arr1" = arr1_GR, "arr10" = arr10_GR, 
                       "arr12" = arr12_GR, "Shoots" = Shoots_all_FC1.2)
vennplot(overlaps_shoots, by = "gplots")

#Plotting gene expression vs. chromatin accessibility CK-induced changes
#Extract geneID and log2FC columns from diffReps
Roots_extracted = Roots_ATAC_genes_all[ , c("geneId", "log2FC")]
Shoots_extracted = Shoots_ATAC_genes_all[ , c("geneId", "log2FC")]

#Extract RNA change for overlapping geneIDs
Roots_extracted$RNA_log2FC = RootsGenes$logFC[match(Roots_extracted$geneId, RootsGenes$gene)]
Shoots_extracted$RNA_log2FC = ShootsGenes$logFC[match(Shoots_extracted$geneId, ShootsGenes$gene)]

#Linear regression
lm = round(summary(lm(formula = Roots_extracted$log2FC ~ Roots_extracted$RNA_log2FC))$r.squared, 2)
lm = round(summary(lm(formula = Shoots_extracted$log2FC ~ Shoots_extracted$RNA_log2FC))$r.squared, 2)

#Plot results
ggplot(Roots_extracted, aes(x = log2FC, y = RNA_log2FC)) + 
  geom_point(colour = "red") + 
  theme_classic() + 
  geom_smooth(method = "lm", colour = "black", size = 1, 
              se = FALSE) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_blank()) + 
  annotate("text", label = paste("R^2: ", 
                                 lm, sep=""), 
           x=0.75, y = -1.5, parse = T, size = 6)

ggplot(Shoots_extracted, aes(x = log2FC, y = RNA_log2FC)) + 
  geom_point(colour = "red") + 
  theme_classic() + 
  geom_smooth(method = "lm", colour = "black", size = 1, 
              se = FALSE) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_blank()) + 
  annotate("text", label = paste("R^2: ", 
                                 lm, sep=""), 
           x=1.0, y = -0.5, parse = T, size = 6)  


