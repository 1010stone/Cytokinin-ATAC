#!/bin/bash -l

# Script to call differential regions

module load perl

diffReps.pl \
  -tr BED_FILES/(Treatment_ATAC_Rep1)_SHIFTED_SORTED.bed \ 
    BED_FILES/(Treatment_ATAC_Rep2).bed \ 
  -co BED_FILES/(Control_ATAC_Rep1)_SHIFTED_SORTED.bed 
    BED_FILES/(ControL_ATAC_Rep2)_SHIFTED_SORTED.bed \
  -ch genome.txt \
  -ns b \
  -mo n \
  -re DIFFREPS/diffReps_WT_Results.txt \
  -me nb"
