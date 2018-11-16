#Load database
tair10 <- loadDb("TAIR10_GFF3_genes_chr_removed.sqlite")
                                                               
#Read in the diffReps info
diffRepsTable <- read.csv("diffRepsTable.csv", stringsAsFactors = F)
diffRepsList <- split(diffRepsTable, f = diffRepsTable$Tissue.Type)
x <- diffRepsTable$Tissue.Type
#Read in diffReps peak files
diffRepsGRList <- lapply(diffRepsList, function(x){
    PEAKS <- readPeakFile(x$Peakfile)
    #-----
    # Name columns to reflect diffreps columns
    names(PEAKS@elementMetadata@listData) <- c("Length", "Treatment Count", "Control Count", "Treatment Ave", "Control Ave", "Treatment Enr", "Control Enr", "Event", "log2FC", "pval", "padj", "winSta", "winEnd", "winFC", "winP", "winQ")
    # Make new columns for metadata about sample
    mcols(PEAKS)$FC <- 2^PEAKS$log2FC
    
    return(PEAKS)
})

#Create annotation of diff regions
diffRepsAnnoList <- lapply(diffRepsGRList, annotatePeak, TxDb=tair10,
  tssRegion = c(-3000, 3000), verbose = TRUE)

#Return all information in a dataframe
diffRepsAnnoList_Dataframe = lapply(diffRepsAnnoList, function(x) as.data.frame(x))

#Create individual data frames
Roots = diffRepsAnnoList_Dataframe$`Roots`
Shoots = diffRepsAnnoList_Dataframe$`Shoots`

#Change first column name
colnames(Roots)[1] = "chr"
colnames(Shoots)[1] = "chr"

#Apply adjusted p-value cutoff
Roots_padj = Roots[Roots$padj <= 0.05, ]
Shoots_padj = Shoots[Shoots$padj <= 0.05, ]

#Apply fold-change cutoff
Roots_FC1.2 = Roots_padj[Roots_padj$FC <= 0.8 | Roots_padj$FC >= 1.2, ]
Shoots_FC1.2 = Shoots_padj[Shoots_padj$FC <= 0.8 | Shoots_padj$FC >= 1.2, ]

#Filter out blacklist region
Roots_FC1.2 = Roots_FC1.2[!(Roots_FC1.2$chr == 2 & Roots_FC1.2$start %in% c(3238000:3629000)), ]
Shoots_FC1.2 = Shoots_FC1.2[!(Shoots_FC1.2$chr == 2 & Shoots_FC1.2$start %in% c(3238000:3629000)), ]

#Separate increased and decreased accessibility regions
Roots_FC1.2_up = Roots_FC1.2[Roots_FC1.2$FC >= 1.2, ]
Roots_FC1.2_down = Roots_FC1.2[Roots_FC1.2$FC <= 0.8, ]
Shoots_FC1.2_up = Shoots_FC1.2[Shoots_FC1.2$FC >= 1.2, ]
Shoots_FC1.2_down = Shoots_FC1.2[Shoots_FC1.2$FC <= 0.8, ]
