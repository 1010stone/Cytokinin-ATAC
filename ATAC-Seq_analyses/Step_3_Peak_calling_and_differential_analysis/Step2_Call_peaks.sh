#!/bin/bash -l

# Script to create tag directory, call peaks, convert peaks to bed file and sort bed file.

module load homer bedtools

ATAC="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    makeTagDirectory \
      $ATAC # output tag directories \
      $ATAC.bed -format bed \ # input is fully processed BED file containing (nuclear, non-duplicated, shifted) reads.
    
    && findPeaks \
      $ATAC \ # output tag directories
      -region -size 200 \ 
      -minDist 202 \
      -o $ATAC_peaks.txt \ # output peak file
      -gsize 135000000 \
      -tbp 0 \
    
    # awk command to process peak.txt file and output genomic coordinates in bed format
    # skips first 35 lines, then outputs genomic coordinates of peaks
    && awk 'NR < 35; NR > 35 {OFS = "\t"; print $2, $3, $4}' $ATAC_peaks.txt \
    # pipes resulting bed file for sorting via bedtools.
    |bedtools sort -i - > $ATAC_peaks.bed  
  done
