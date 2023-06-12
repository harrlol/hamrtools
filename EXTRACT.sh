#!/bin/bash
set -eu

# Suite for Modification or Annotation-purposed Cleaning and Keeping

if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: EXTRACT.sh <proj dir> <knownant> <distTECH> <distGENO> <gene annotation file>"
echo "EXAMPLE:"
exit 1
fi

# Project directory, the same one you fed into SMACK
dir=$1

# The csv (in modtbl format) of the known mod you want analyzed in distToKnownMod
antcsv=$2

# The seq_tech of the sample group you want analyzed in distToKnownMod
tech=$3

# The genotype of the sample group you want analyzed in distToKnownMod
geno=$4

# The gene annotation file for your organism, obtained with Diep's code
ant=$5
curdir=$(dirname $0)

# collapse all overlapped data into longdf
Rscript $curdir/allLapPrep.R \
    $dir/lap \
    $dir

# overview of modification proportaion
Rscript $curdir/abundByLap.R \
    $dir/mod_long.csv \
    gene \
    $dir

# analyze hamr-mediated/true clustering across project
Rscript $curdir/clusterAnalysis.R \
    $dir/mod_long.csv \
    $dir

# analyze hamr-mediated/true clustering across project
Rscript $curdir/distToKnownMod.R \
    $dir/mod_long.csv \
    $antcsv \
    gene $tech $geno \
    $dir

# looking at around which part of a gene a certain mod is likely to be found 
Rscript $curdir/modDistribution.R \
    $dir/mod_long.csv \
    $ant \
    $dir