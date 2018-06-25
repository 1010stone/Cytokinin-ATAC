#!/bin/bash -l

# Script to subsample 20 million reads from each BAM file and index them.

# Load in modules
module load samtools bbmap bedtools

# Create variable containing file names.
ATAC=(ATAC_1 ATAC_2 ATAC_n)

# Create variable containing number of reads to subsample. 
# This script was also used to subsample 5 M reads to obtain SPOT scores.
READS="20000000"

# Create variable that contains the total number of reads for each dataset.
declare -a total=()

for i in "${!ATAC[@]}"
do
	total+=( "$(samtools view -c -F 260 BAM_FILES/${ATAC[i]}_SHIFTED_SORTED.bam)" )
	echo ${arr[*]}
done

declare -p total


#Calculates the normalization factor to get to 20M reads per library.
declare -a factor=()


for i in "${!total[@]}"
do
	factor+=($(bc <<< "scale=2; $READS/${total[i]}"))
	echo ${factor[*]}
done

declare -p factor


for i in ${!ATAC[@]}
do
	samtools view -s ${factor[i]} -b BAM_FILES/${ATAC[i]}_SHIFTED_SORTED.bam >BAM_FILES/${ATAC[i]}_SUBSAMPLED.bam" \
	&& samtools index BAM_FILES/${ATAC[i]}_SUBSAMPLED.bam
done
