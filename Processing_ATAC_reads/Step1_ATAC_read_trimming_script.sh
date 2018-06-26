#!/bin/bash -l

# Example bash script used to trim reads in each of the ATAC files

# Load in the bbmap module
module load bbmap

# Determine the bbmap version number \
# this is needed to find the correct directory \
# containing the updated Arabidopsis index \
# for correct version of bbmap (it gets updated frequently...)
bbmap=`echo $LOADEDMODULES| awk -F'bbmap/' '{print \$2}'`

# Create variable for unique directory names that correspond to different library names.
NAMES=(BAM_1 BAM_2 BAM_n)

ATAC=(ATAC_1 ATAC_2 ATAC_n)

for i in ${!ATAC[@]}
do
  bbduk.sh \
    in=${ATAC[i]}_R1_001.fastq.gz \ # input read pair 1
    in2=${ATAC[i]}_R2_001.fastq.gz \ # input read pair 2
    out=${NAMES[i]}-R1.trimmed.fastq \ # output trimmed reads for read pair 1
    out2=${NAMES[i]}-R2.trimmed.fastq \ # output trimmed reads for read pair 2
    ref=/nas02/apps/bbmap-$bbmap/bbmap/resources/adapters.fa \ # file containing Illumina adapter sequences
    k=13 \ # kmers to trim
    ktrim=r # direction of ktrimming
done
