#Load in packages
library(agricolae)
library(dplyr)
library(lsmeans)
library(edgeR)
library(DESeq2)
library(multcomp)
library(multcompView)
library(DTK)
library(DescTools)
library(dunn.test)
library(gdata)
library(xlsx)
library(readxl)
library(Rsubread)
library(ggplot2)
library(gplots)
library(multtest)
library(VennDiagram)
library(reshape2)
library(RUVSeq)

# Roots RUV-Seq 
#Differential Gene Expression Analysis
list = list("COUNTS/1_ROOTS_NAOH_COUNTS.txt", "COUNTS/2_ROOTS_NAOH_COUNTS.txt", 
            "COUNTS/3_ROOTS_NAOH_COUNTS.txt", "COUNTS/1_ROOTS_BA_COUNTS.txt", 
            "COUNTS/2_ROOTS_BA_COUNTS.txt", "COUNTS/3_ROOTS_BA_COUNTS.txt")
groups = c("Roots_NaOH", "Roots_NaOH", "Roots_NaOH", "Roots_BA", "Roots_BA", "Roots_BA")
labels = c("Roots_NaOH_1", "Roots_NaOH_2", "Roots_NaOH_3", "Roots_BA_1", "Roots_BA_2", "Roots_BA_3")
counts = readDGE(list, path=NULL, columns=c(1,7), group = groups, labels = labels, header = TRUE, skip = 1)
# Subset the root data
group <- factor(rep(c("Roots-Control", "Roots-BA"),each=3),levels=c("Roots-Control", "Roots-BA"))
group

# Filter data by CPM using edgeR
library(edgeR)
raw <- DGEList(counts=counts,group=group)
nrow(raw)
keep <- rowSums(cpm(raw)>1) >= 3 # at least 1 CPM in 3 or more samples 
filtered <- raw[keep, keep.lib.sizes=T]
nrow(filtered) #19117

###   RUVseq  #####
##################
library(RUVSeq)
genes <- rownames(filtered)
colnames(filtered)
library(RColorBrewer)

# Before RUV
## Before normalization
set <- newSeqExpressionSet(as.matrix(filtered), 
                           phenoData = data.frame(group, row.names =colnames(filtered)))
set

## After normalization
set1 <- betweenLaneNormalization(set, which="upper")
colors<-c("#40004b","#762a83","#e7d4e8",
          "#9970ab","#9970ab","#9970ab")

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[group])
par(bg = "white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9, 11), xpd=TRUE,ps =13)
plotPCA(set1, col=colors, 
        cex=2, labels=F, pch = 17)
legend("topright",cex=1.2,pch = 17, 
       legend = c("Roots-Control-1", "Roots-Control-2", "Roots-Control-3", 
                  "Roots-BA-1", "Roots-BA-2", "Roots-BA-3"), 
       col=colors,inset=c(-0.37,0),bty = "n")

#RUVs: Estimating the factors of unwanted variation using replicate samples
differences<-makeGroups(group)
differences

set3 <- RUVs(set1, genes, k=1, differences)
pData(set3)

colors <- brewer.pal(4, "Set1")

plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[group])
par(bg = "white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9, 11), xpd=TRUE,ps =13)
plotPCA(set3, col=colors[group], 
        cex=2, labels=F, pch = 17,
        main="RUVs normalization - replicate samples")
legend("topright",cex=1.2,pch = 17, 
       legend = c("Roots-BA","Roots"), 
       col=colors,inset=c(-0.25,0),bty = "n")

### edgeR ####
##############
design <- model.matrix(~group+W_1,data=pData(set3))
design

colnames(design)
y <- DGEList(counts=counts(set3),group=group)
y <- calcNormFactors(y, method = 'upper')

y$samples
levels(y$samples$group)

#plotMDS(y)

y<-estimateDisp(y,design,robust = T)
y$common.dispersion #0.0135056

plotBCV(y)

fit <-glmFit(y, design)

########################
#### WT BA response ####
########################
lrt <-glmLRT(fit, coef=2) # compare WT-NaOH and WT-BA
top <- topTags(lrt, n=nrow(set3))$table
sum(top$FDR<0.05 & abs(top$logFC)>0.584936) #979

RootsGenes = top[top$FDR<0.05 & abs(top$logFC)>0.584936, ]
RootsGenes$geneId = row.names(RootsGenes)
row.names(RootsGenes) = c()

RootsGenes_up = RootsGenes[RootsGenes$logFC > 0.5, ]
RootsGenes_down = RootsGenes[RootsGenes$logFC < 0.5, ]






#Shoots RUV-Seq
#Differential Gene Expression Analysis
list = list("COUNTS/1_SHOOTS_NAOH_COUNTS.txt", "COUNTS/2_SHOOTS_NAOH_COUNTS.txt", 
            "COUNTS/3_SHOOTS_NAOH_COUNTS.txt", "COUNTS/1_SHOOTS_BA_COUNTS.txt", 
            "COUNTS/2_SHOOTS_BA_COUNTS.txt", "COUNTS/3_SHOOTS_BA_COUNTS.txt")
groups = c("SHOOTS_NaOH", "SHOOTS_NaOH", "SHOOTS_NaOH", "SHOOTS_BA", "SHOOTS_BA", "SHOOTS_BA")
labels = c("SHOOTS_NaOH_1", "SHOOTS_NaOH_2", "SHOOTS_NaOH_3", "SHOOTS_BA_1", "SHOOTS_BA_2", "SHOOTS_BA_3")
counts = readDGE(list, path=NULL, columns=c(1,7), group = groups, labels = labels, header = TRUE, skip = 1)
# Subset the shoot data
group <- factor(rep(c("SHOOTS-Control", "SHOOTS-BA"),each=3),levels=c("SHOOTS-Control", "SHOOTS-BA"))
group

# Filter data by CPM using edgeR
library(edgeR)
raw <- DGEList(counts=counts,group=group)
nrow(raw)
keep <- rowSums(cpm(raw)>1) >= 3 # at least 1 CPM in 3 or more samples 
filtered <- raw[keep, keep.lib.sizes=T]
nrow(filtered) #18428

###   RUVseq  #####
##################
library(RUVSeq)
genes <- rownames(filtered)
colnames(filtered)
library(RColorBrewer)

# Before RUV
## Before normalization
set <- newSeqExpressionSet(as.matrix(filtered), 
                           phenoData = data.frame(group, row.names =colnames(filtered)))
set

## After normalization
set1 <- betweenLaneNormalization(set, which="upper")
colors<-c("#40004b","#762a83","#e7d4e8",
          "#9970ab","#9970ab","#9970ab")

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[group])
par(bg = "white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9, 11), xpd=TRUE,ps =13)
plotPCA(set1, col=colors, 
        cex=2, labels=F, pch = 17)
legend("topright",cex=1.2,pch = 17, 
       legend = c("SHOOTS-Control-1", "SHOOTS-Control-2", "SHOOTS-Control-3", 
                  "SHOOTS-BA-1", "SHOOTS-BA-2", "SHOOTS-BA-3"), 
       col=colors,inset=c(-0.37,0),bty = "n")

#RUVs: Estimating the factors of unwanted variation using replicate samples
differences<-makeGroups(group)
differences

set3 <- RUVs(set1, genes, k=1, differences)
pData(set3)

colors <- brewer.pal(4, "Set1")

plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[group])
par(bg = "white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9, 11), xpd=TRUE,ps =13)
plotPCA(set3, col=colors[group], 
        cex=2, labels=F, pch = 17,
        main="RUVs normalization - replicate samples")
legend("topright",cex=1.2,pch = 17, 
       legend = c("SHOOTS-BA","SHOOTS"), 
       col=colors,inset=c(-0.25,0),bty = "n")

### edgeR ####
##############
design <- model.matrix(~group+W_1,data=pData(set3))
design

colnames(design)
y <- DGEList(counts=counts(set3),group=group)
y <- calcNormFactors(y, method = 'upper')

y$samples
levels(y$samples$group)

#plotMDS(y)

y<-estimateDisp(y,design,robust = T)
y$common.dispersion #0.01326261

plotBCV(y)

fit <-glmFit(y, design)

########################
#### WT BA response ####
########################
lrt <-glmLRT(fit, coef=2) # compare WT-NaOH and WT-BA
top <- topTags(lrt, n=nrow(set3))$table
sum(top$FDR<0.05 & abs(top$logFC)>0.584936) #337

ShootsGenes = top[top$FDR<0.05 & abs(top$logFC)>0.584936, ]
ShootsGenes$geneId = row.names(ShootsGenes)
row.names(ShootsGenes) = c()

ShootsGenes_up = ShootsGenes[ShootsGenes$logFC > 0.5, ]
ShootsGenes_down = ShootsGenes[ShootsGenes$logFC < 0.5, ]
