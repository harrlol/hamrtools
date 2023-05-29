#!/bin/bash
set -eu

# Script to use sratool (newest version) that is otherwise unavailable in bash 

# First, install the relative sratoolkit manually from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

# Next, create an SRR accession list (in .txt format) using the SRA run selector tool.

# Locate both your sratoolkit and your .txt


if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: fasterq-dumpAdaptor.sh <sratoolkit dir> <accession.txt> <out dir>"
echo "EXAMPLE:"
exit 1
fi

dir=$1
txt=$2
out=$3

if [ ! -d "$out/cleaned" ] 
then
    echo "You can find your cleaned fastq files at $out/cleaned" 
    mkdir $out/cleaned
fi

if [ ! -d "$out/trimmed" ] 
then
    echo "You can find your trimmed fastq files at $out/trimmed" 
    mkdir $out/trimmed
fi

if [ ! -d "$out/fastqc_results" ] 
then
    echo "You can find all the fastqc test results at $out/fastqc_results" 
    mkdir $out/fastqc_results
fi

while IFS= read -r line
do
    echo "begin downloading $line..."
    $dir/bin/fasterq-dump "$line" -O $out --verbose
    echo "combining $line into a single fastq file..."
    cat $out/$line*.fastq > $out/cleaned/$line.fastq
    echo "$line complete, trimming..."
    echo ""
    ~/TrimGalore-0.6.10/trim_galore -o $out/trimmed $out/cleaned/$line.fastq
    echo "trimming complete, performing fastqc..."
    echo ""
    fastqc $out/trimmed/$line"_trimmed.fq" -o $out/fastqc_results &
    fastqc $out/cleaned/$line.fastq -o $out/fastqc_results
    echo "$line complete"
    echo ""
done < "$txt"