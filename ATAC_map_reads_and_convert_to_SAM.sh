#!/bin/bash -l

# Example script used to convert SAM to BAM and sort and index BAM files

# Load in modules
module load samtools bbmap bowtie deeptools bedtools

# Create variable that contains the unique IDs for names of directories/files.
NAMES="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    bowtie -X 1000 -t -p 4 -v 2 --best --strata -m 2 -y -k 1 \
      /proj/seq/data/TAIR10_Ensembl/Sequence/BowtieIndex/genome \ 
      -1 FASTQ_FILES/$ATAC-R1.trimmed.fastq \ 
      -2 FASTQ_FILES/$ATAC-R2.trimmed.fastq \ 
      -S SAM_FILES/$ATAC_MAPPED.sam \
    && samtools view -b SAM_FILES/$ATAC_MAPPED.sam \ 
      |samtools sort - -o BAM_FILES/$ATAC_MAPPED_SORTED.bam \
    && samtools index BAM_FILES/$POTTER_MAPPED_SORTED.bam
  done
