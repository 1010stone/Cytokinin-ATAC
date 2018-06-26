#!/bin/bash -l

# Script to count reads within differential regions and plot averaged read densities.

module load bedtools

bedtools multicov \
  -bams ATAC_merged_treated_datasets.bam \
        ATAC_merged_control_datasets.bam \
  -bed Differential_regions.bed \
  >matrix.txt
