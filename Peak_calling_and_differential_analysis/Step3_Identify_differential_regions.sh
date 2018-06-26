#!/bin/bash -l

# Script to call differential regions

module load perl

diffReps.pl \
  -tr ATAC_Treatment_Rep1.bed \ 
      ATAC_Treatment_Rep2.bed \
      ATAC_Treatment_RepN.bed \
  -co ATAC_Control_Rep1.bed \ 
      ATAC_Control_Rep2.bed \
      ATAC_Control_RepN.bed
  -ch genome.txt \
  -ns b \
  -mo n \
  -re output_results.txt \
  -me nb"
