#!/bin/bash
set -eu

# Script to overlap consensus mods with various libraries (UTR, CDS, gene, etc.) after post HAMR processing in R.
# the libraries are obtained in various banks like TAIR

if [ "$#" -lt 6 ]; then
echo "Missing arguments!"
echo "USAGE: consensusOverlap.sh <consensus.bed> <cds.bed> <utr.bed> <gene.bed> <mrna.bed> <out dir>"
echo "EXAMPLE:"
exit 1
fi

smp=$1
cds=$2
utr=$3
gene=$4
mrna=$5
out=$6


IFS="/" read -ra sections <<< "$smp"
temp="${sections[-1]}"

IFS="." read -ra templ <<< "$temp"
smpname="${templ[0]}"

echo ""
echo "sample = $1"
echo "CDS lib = $2"
echo "UTR lib = $3"
echo "Gene lib = $4"
echo "mRNA lib = $5"
echo "out = $6"
echo "consensus file prefix: $smpname"
echo ""

echo "..."
#overlap with cds
intersectBed \
    -a $cds \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_CDS".bed
echo "finished finding overlap with CDS library"

echo "..."
#overlap with utr
intersectBed \
    -a $utr \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_UTR".bed
echo "finished finding overlap with UTR library"

echo "..."
#overlap with gene
intersectBed \
    -a $gene \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_gene".bed
echo "finished finding overlap with gene library"

echo "..."
#overlap with mrna
intersectBed \
    -a $mrna \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_mRNA".bed
echo "finished finding overlap with mRNA library"