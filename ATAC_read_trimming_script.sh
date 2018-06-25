#!/bin/bash -l

# Example bash script used to trim reads in each of the ATAC files

# Load in the bbmap module
module load bbmap

# Determine the bbmap version number \
# this is needed to find the correct directory \
# containing the updated Arabidopsis index
bbmap=`echo $LOADEDMODULES| awk -F'bbmap/' '{print \$2}'`

# Create variable for unique directory names that correspond to different library names.
NAMES=(BAM_1 BAM_2 BAM_n)

ATAC=(ATAC_1_CGAGGCTG_S1_L008 ATAC_2_AAGAGGCA_S2_L008 ATAC_3_TAAGGCGA_S3_L008 ATAC_4_CGTACTAG_S4_L008)

for i in ${!ATAC[@]}
do
  bbduk.sh in=../../../HTSF/161202_UNC32-K00270_0031_BHCKF2BBXX/${ATAC[i]}_R1_001.fastq.gz \
    in2=../../../HTSF/161202_UNC32-K00270_0031_BHCKF2BBXX/${ATAC[i]}_R2_001.fastq.gz \
    out=FASTQ_FILES/${NAMES[i]}-R1.trimmed.fastq \
    out2=FASTQ_FILES/${NAMES[i]}-R2.trimmed.fastq \
    ref=/nas02/apps/bbmap-$bbmap/bbmap/resources/adapters.fa \
    k=13 \
    ktrim=r
done
