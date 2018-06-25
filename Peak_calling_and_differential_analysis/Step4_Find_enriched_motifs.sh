#!/bin/bash -l

# Script to find enriched motifs in differential regions.

module load homer

findMotifsGenome.pl \
  DIFFREPS/(Differtial_regions).bed \
  tair10 \
  DIFFREPS/MOTIFS/ \
  -size given \
  -len 6,10,12 \
  -preparse \
