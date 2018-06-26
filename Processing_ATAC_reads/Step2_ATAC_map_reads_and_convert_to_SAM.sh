#!/bin/bash -l

# Example script used to map reads, convert SAM to BAM and then sort and index BAM files

# Load in modules
module load samtools bbmap bowtie deeptools bedtools

# Create variable that contains the unique IDs for names of directories/files.
NAMES="ATAC_1 ATAC_2 ATAC_n"

for ATAC in $ATAC
  do
    bowtie -X 1000 -t -p 4 -v 2 --best --strata -m 2 -y -k 1 \
      /proj/seq/data/TAIR10_Ensembl/Sequence/BowtieIndex/genome \ # location of bowtie index for genome
      -1 $ATAC-R1.trimmed.fastq \ # input trimmed reads for read pair 1
      -2 $ATAC-R2.trimmed.fastq \ # input trimmed reads for read pair 2
      -S $ATAC_MAPPED.sam \ # output mapped reads in sam format
   
    # convert sam to bam
    && samtools view -b $ATAC_MAPPED.sam \ 
      |samtools sort - -o $ATAC_MAPPED_SORTED.bam \
    
    # index bam file
    && samtools index $ATAC_MAPPED_SORTED.bam
  done
