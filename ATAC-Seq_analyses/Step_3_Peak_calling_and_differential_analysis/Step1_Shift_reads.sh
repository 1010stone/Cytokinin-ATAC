#!/bin/bash -l

# Script to shift reads +4 on forward and -5.

# Load in modules
module load samtools bedtools

# Create variable containing unique descriptors of your datasets.
ATAC="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    awk 'BEGIN {OFS ="\t"} {if ($6 == "+"){print $1, $2 + 4, $3 + 4, $4, $5, $6} \
      else{print $1, $2 - 5, $3 - 5, $4, $5, $6}}' \
      $ATAC.bed >$ATAC_SHIFTED_UNSORTED.bed \
    
    # now sort resulting bed file
    && sort -k1,1 -k2,2n -k3,3n $ATAC_SHIFTED_UNSORTED.bed \
      >$ATAC_SHIFTED_SORTED.bed \
    
    # convert bed file back to bam file
    && bedtools bedtobam -i $ATAC_SHIFTED_SORTED.bed \
      -g genome.txt >$ATAC_SHIFTED_SORTED.bam \
    
    # index new bam containing shifted and sorted reads
    && samtools index $ATAC_SHIFTED_SORTED.bam
  done
