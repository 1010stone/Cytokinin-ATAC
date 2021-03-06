#!/bin/bash -l

# Script to find enriched motifs in differential regions.

module load homer

findMotifsGenome.pl \
  Differtial_regions.bed \
  tair10 \
  output_directory/ \
  -size given \
  -len 6,10,12 \
  -preparse \
