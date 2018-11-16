#!/bin/bash -l
# tophat_manyfiles.sbatch
#
#SBATCH -p general
#SBATCH -t 2-00:00
#SBATCH --mem=24000

module load samtools tophat subread

now=$(date +"%m_%d_%Y")

echo $now
echo $1 $2 $3
# script to map, index and count reads
genes="/proj/seq/data/MSU6_rice_Ensembl/Annotation/Genes/genes.gtf"
ref="/proj/seq/data/MSU6_rice_Ensembl/Sequence/Bowtie2Index/genome"

echo "mapping ${3}"
tophat -p 4 -G $genes -o Results/$3.thout $ref $1 $2_2.fq.gz

echo "indexing ${3}"
samtools index Results/$3.thout/accepted_hits.bam

echo "counting feature counts for ${3}"
featureCounts -t exon -g gene_id -a $genes -o Counts/$3_Counts.txt Results/$3.thout/accepted_hits.bam

echo "All done :-)"

