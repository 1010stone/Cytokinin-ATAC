#!/bin/bash -l

# Script to merged biological replicates.

samtools merge \
  ATAC_(Control or Treated)_MERGED.bam \ # output of merged bam file
  ATAC_1_SHIFTED_SORTED.bam \ # input bam 1
  ATAC_2_SHIFTED_SORTED.bam \ # input bam 2
  ATAC_n_SHIFTED_SORTED.bam \ # input bam n
