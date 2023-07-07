#!/bin/bash
set -u

# Script to use sratool (newest version) that is otherwise unavailable in bash 

# First, install the relative sratoolkit manually from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

# Next, create an SRR accession list (in .txt format) using the SRA run selector tool.

# Locate both your sratoolkit and your .txt


if [ "$#" -lt 3 ]; then
echo "Missing arguments!"
echo "USAGE: fasterq-dumpAdaptor.sh <out dir> <T/F> <SRR>"
echo "EXAMPLE:"
exit 1
fi

line=$1
out=$2
trim=$3

echo "begin downloading $line..."
fasterq-dump "$line" -O $out/raw --verbose

echo "[$line] performing fastqc on raw file..."
fastqc $out/raw/$line.fastq -o $out/fastqc_results &

if [ ! "$trim" = 'T' ]; then
    echo "[$line] trimming..."
    trim_galore -o $out/trimmed $out/raw/$line.fastq

    echo "[$line] trimming complete, performing fastqc..."
    fastqc $out/trimmed/$line"_trimmed.fq" -o $out/fastqc_results
fi

echo "finished processing $line"
echo ""