#!/bin/bash -l

# Script to compare correlation between BAM files.

module load deeptools

multiBamSummary bins --bamfiles \
  ATAC_1_SUBSAMPLED.bam \
  ATAC_2_SUBSAMPLED.bam \
  ATAC_n_SUBSAMPLED.bam \
  -out results.npz \
  --outRawCounts readCounts.tab \

&& plotCorrelation --corData results.npz \
--corMethod spearman --whatToPlot heatmap --skipZeros --removeOutliers \
--plotFile heatmap.png --plotNumbers \
--outFileCorMatrix HeatmapSpearmanCorr_readCounts.tab \
--colorMap Reds \
--labels ATAC_1 ATAC_2 ATAC_n \
--zMin 0.70 \
--zMax 1.00
