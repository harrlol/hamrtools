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
fiveutr=$3
threeutr=$4
gene=$5
mrna=$6
exon=$7
all=$8
out=$9

IFS="/" read -ra sections <<< "$smp"
temp="${sections[-1]}"

IFS="." read -ra templ <<< "$temp"
smpname="${templ[0]}"

echo ""
echo "sample = $1"
echo "out = $9"
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
#overlap with 5utr
intersectBed \
    -a ${fiveutr} \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_fiveUTR".bed
echo "finished finding overlap with 5UTR library"

echo "..."
#overlap with 3utr
intersectBed \
    -a ${threeutr} \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_threeUTR".bed
echo "finished finding overlap with 3UTR library"

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
    > $out/$smpname"_primarymRNA".bed
echo "finished finding overlap with primary mRNA library"

echo "..."
#overlap with exon
intersectBed \
    -a $exon \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_exon".bed
echo "finished finding overlap with exon library"

echo "..."
#overlap with nc rna
intersectBed \
    -a $all \
    -b $smp \
    -wa -wb \
    > $out/$smpname"_ncRNA".bed
echo "finished finding overlap with ncRNA library"