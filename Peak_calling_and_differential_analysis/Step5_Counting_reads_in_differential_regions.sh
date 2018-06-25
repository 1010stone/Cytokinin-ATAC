#!/bin/bash -l

# Script to count reads within differential regions.

module load samtools bedtools

bedtools multicov \
  -bams \
    BAM_FILES/ATAC_1.bam \
    BAM_FILES/ATAC_2.bam \
    BAM_FILES/ATAC_n.bam \
  -bed DIFFREPS/Roots_WT_All_diffRepsFC1.2.bed \
  >DIFFREPS/ATAC_COMPARING_ROOTS_REGIONS_FILE.txt
