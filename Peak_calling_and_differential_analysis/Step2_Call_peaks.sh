#!/bin/bash -l

# Script to create tag directory, call peaks, convert peaks to bed file and sort bed file.

module load homer bedtools

ATAC="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    makeTagDirectory TAG_DIRECTORIES/$ATAC BED_FILES/$ATAC_SHIFTED_SORTED.bed -format bed \
    
    && findPeaks TAG_DIRECTORIES/$ATAC -region -size 200 \ 
    -minDist 202 -o PEAK_FILES/$ATAC_peaks.txt -gsize 135000000 -tbp 0 \
    
    && awk 'NR < 35; NR > 35 {OFS = "\t"; print $2, $3, $4}' PEAK_FILES/$ATAC_peaks.txt \ 
    |bedtools sort -i - >PEAK_FILES/$ATAC_peaks.bed  
  done
