#!/bin/bash
set -eu

# Simply extract out information regarding HAMR predictions

if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: EXTRACT.sh <proj dir> <knownant> <distTECH> <distGENO> <gene annotation file>"
echo "EXAMPLE:"
exit 1
fi

# Project directory, the same one you fed into SMACK
dir=$1

# The csv (in modtbl format) of the known mod you want analyzed in distToKnownMod
antcsv=$2

# The annotation folder for your organism
ant=$3
curdir=$(dirname $0)

# collapse all overlapped data into longdf
Rscript $curdir/allLapPrep.R \
    $dir/lap \
    $dir

# overview of modification proportion
Rscript $curdir/abundByLap.R \
    $dir/mod_long.csv \
    $ant \
    $dir

# analyze hamr-mediated/true clustering across project
Rscript $curdir/clusterAnalysis.R \
    $dir/mod_long.csv \
    $dir

# analyze hamr-mediated/true clustering across project
Rscript $curdir/distToKnownMod.R \
    $dir/mod_long.csv \
    $antcsv

gff=$(find $ant -maxdepth 1 -name "*_gene*")
if [ -z "$gff" ]; then
    echo "gff3 annotation file not found"
    exit 1
fi

# looking at around which part of a gene a certain mod is likely to be found
Rscript $curdir/modDistribution.R \
    $dir/mod_long.csv \
    $gff \
    $dir

# looking at RNA subtype for mods
Rscript $curdir/RNAtype.R \
    $dir/mod_long.csv