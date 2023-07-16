#!/bin/bash
set -eu

# Simply extracts out information regarding HAMR predictions

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: EXTRACT.sh <proj dir> <gene annotation> [knownant] "
echo "EXAMPLE:"
exit 1
fi

# Project directory, the same one you fed into SMACK
dir=$1

# The annotation folder for your organism
ant=$2

curdir=$(dirname $0)

panther=$3

echo "activating rspace..."
eval "$(conda shell.bash hook)"
conda activate rspace
wait
echo "rspace activated"
echo ""

echo "generating long modification table..."
# collapse all overlapped data into longdf
Rscript $curdir/allLapPrep.R \
    $dir/lap \
    $dir
echo "done"
echo ""

echo "plotting modification abundance..."
# overview of modification proportion
Rscript $curdir/abundByLap.R \
    $dir/mod_long.csv \
    $ant \
    $dir
echo "done"
echo ""

echo "performing modification cluster analysis..."
# analyze hamr-mediated/true clustering across project
Rscript $curdir/clusterAnalysis.R \
    $dir/mod_long.csv \
    $dir
echo "done"
echo ""

if [ ! -z "${4+x}" ]; then
    echo "known modification landscape provided, performing relative positional analysis to known mod..."
    # The csv (in modtbl format) of the known mod you want analyzed in distToKnownMod
    antcsv=$4
    # analyze hamr-mediated/true clustering across project
    Rscript $curdir/distToKnownMod.R \
        $dir/mod_long.csv \
        $antcsv
    echo "done"
    echo ""
else 
    echo "known modification file not detected, skipping relative positional analysis"
    echo ""
fi

gff=$(find $ant -maxdepth 1 -name "*_gene*")
if [ -z "$gff" ]; then
    echo "gff3 annotation file not found"
    exit 1
fi

echo "classifying modified RNA subtype..."
# looking at RNA subtype for mods
Rscript $curdir/RNAtype.R \
    $dir/mod_long.csv
echo "done"
echo ""

if [ ! -d "$dir/go" ]; then mkdir $dir/go; echo "created path: $dir/go"; fi

if [ ! -d "$dir/go/genelists" ]; then mkdir $dir/go/genelists; echo "created path: $dir/go/genelists"; fi

if [ ! -d "$dir/go/pantherout" ]; then mkdir $dir/go/pantherout; echo "created path: $dir/go/pantherout"; fi

echo "performing gene ontology analysis..."
$curdir/modGO.sh \
    $dir \
    $panther
echo "done"

echo "classifying modified RNA subtype..."
# looking at RNA subtype for mods
Rscript $curdir/RNAtype.R \
    $dir/mod_long.csv
echo "done"
echo ""

c=$(find $ant -type f -name "*_CDS.bed")
f=$(find $ant -type f -name "*_fiveUTR.bed")
t=$(find $ant -type f -name "*_threeUTR.bed")
echo "mapping modification regional distribution landscape..."
# looking at RNA subtype for mods
Rscript $curdir/modRegionMapping.R \
    $dir/mod_long.csv \
    $f \
    $c \
    $t
echo "done"
echo ""