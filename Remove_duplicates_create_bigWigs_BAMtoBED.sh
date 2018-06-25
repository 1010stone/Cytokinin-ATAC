#!/bin/bash -l

# Script to remove duplicate reads, convert BAM to BED and make bigWigs from BAM files.

# Load in modules
module load samtools bedtools deeptools

# Create variable containing unique names for directories/file names.
ATAC="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    samtools rmdup BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED.bam \ 
      BAM_FILES/$ATAC"_MAPPED_SORTED_EXTRACTED_CLEANED.bam \
    
    && samtools index BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED_CLEANED.bam \
    
    && bedtools bamtobed -i BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED_CLEANED.bam \
      >BED_FILES/$ATAC.bed \
    
    && samtools idxstats BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED_CLEANED.bam \
      >BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED_CLEANED.idxstats.txt \
    
    && bamCoverage -o BIGWIGS/$ATAC_allReads.bw \ 
      -b BAM_FILES/$ATAC_MAPPED_SORTED_EXTRACTED_CLEANED.bam -of 'bigwig' -p 8 \ 
       --normalizeTo1x 135000000 --binSize 10 \
  done
