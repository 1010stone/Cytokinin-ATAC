#!/bin/bash -l

# Script to merged biological replicates.

samtools merge \
  BAM_FILES/ATAC_(Control or Treated)_MERGED.bam \
  BAM_FILES/ATAC_1_SHIFTED_SORTED.bam \
  BAM_FILES/ATAC_2_SHIFTED_SORTED.bam \
  BAM_FILES/ATAC_n_SHIFTED_SORTED.bam \
