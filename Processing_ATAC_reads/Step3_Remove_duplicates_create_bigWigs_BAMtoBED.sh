#!/bin/bash -l

# Script to remove duplicate reads, convert BAM to BED and make bigWigs from BAM files.

# Load in modules
module load samtools bedtools deeptools

# Create variable containing unique names for directories/file names.
ATAC="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    # remove duplicates
    samtools rmdup $ATAC_MAPPED_SORTED.bam \ 
      $ATAC_MAPPED_SORTED_CLEANED.bam \
    
    # index bam file
    && samtools index $ATAC_MAPPED_SORTED_CLEANED.bam \
    
    # convert bam to bed
    && bedtools bamtobed -i $ATAC_MAPPED_SORTED_CLEANED.bam \
      >$ATAC.bed \
    
    # create bigwig file
    && bamCoverage -o $ATAC_allReads.bw \ 
      -b $ATAC_MAPPED_SORTED_CLEANED.bam -of 'bigwig' -p 8 \ 
       --normalizeTo1x 135000000 --binSize 10 \
  done
