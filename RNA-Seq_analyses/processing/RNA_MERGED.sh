#!/bin/bash -l

# Script to merge sam files and make BigWig files.

names="ROOTS_BA ROOTS_NAOH SHOOTS_BA SHOOTS_NAOH"

module load samtools deeptools

for names in $names
do
samtools merge MERGED/$names.merged.bam file_1.bam file_2.bam file_n.bam \
samtools index MERGED/$names.merged.bam 
bamCoverage -o MERGED/$names.allReads.bw -b MERGED/$names.merged.bam -of 'bigwig' -p 8 --normalizeTo1x 135000000 --binSize 10"
done
