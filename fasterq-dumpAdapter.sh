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


while IFS= read -r line
do
    echo "begin downloading $line"
    $dir/bin/fasterq-dump "$line" -O $out --verbose
    echo "combining $line into a single fastq file"
    cat $out/$line*.fastq > $out/$line.fastq
    echo "$line complete"
    echo ""
done < "$txt"