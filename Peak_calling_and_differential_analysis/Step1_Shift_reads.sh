#!/bin/bash -l

# Script to shift reads +4 on forward and -5.

ATAC="ATAC_1 ATAC_2 ATAC_n"

module load samtools bedtools

for ATAC in $ATAC
  do
    awk 'BEGIN {OFS ="\t"} {if ($6 == "+"){print $1, $2 + 4, $3 + 4, $4, $5, $6} \
      else{print $1, $2 - 5, $3 - 5, $4, $5, $6}}' \
      BED_FILES/$ATAC.bed >BED_FILES/$ATAC_SHIFTED_UNSORTED.bed \
    
    && sort -k1,1 -k2,2n -k3,3n BED_FILES/$ATAC_SHIFTED_UNSORTED.bed \
      >BED_FILES/$ATAC_SHIFTED_SORTED.bed \
    
    && bedtools bedtobam -i BED_FILES/$ATAC_SHIFTED_SORTED.bed \
      -g genome.txt >BAM_FILES/$ATAC_SHIFTED_SORTED.bam \
    
    && samtools index BAM_FILES/$ATAC_SHIFTED_SORTED.bam
  done
