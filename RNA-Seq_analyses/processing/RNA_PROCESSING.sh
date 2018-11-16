#!/bin/bash -l

mkdir Results Counts

now=$(date +"%m_%d_%Y")

echo $now

READ_1=$(ls FASTQ_FILES/*_1.fq.gz)

for f in $READ_1
  do
    READ_2=${f%_1.fq.gz}
    prefix=${f%/*_1.fq.gz}
    echo $f
    echo $READ_2
    echo $prefix
    sbatch -p general --mem=20g -t 02-00:00 tophat.sbatch $f $READ_2 $prefix
done
