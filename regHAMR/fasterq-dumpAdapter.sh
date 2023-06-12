#!/bin/bash
set -eu

# Script to use sratool (newest version) that is otherwise unavailable in bash 

# First, install the relative sratoolkit manually from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

# Next, create an SRR accession list (in .txt format) using the SRA run selector tool.

# Locate both your sratoolkit and your .txt


if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "USAGE: fasterq-dumpAdaptor.sh <sratoolkit dir> <accession.txt> <out dir> <SE/PE> <T/F>"
echo "EXAMPLE:"
exit 1
fi

dir=$1
txt=$2
out=$3
paired_end_flag=$4
pretrimmed_flag=$5


# Create directory to store cleaned (contacenated) fastq files if paired-end not specified
if [ ! "$paired_end_flag" = 'PE' ]; then
    if [ ! -d "$out/cleaned" ] 
    then
        mkdir $out/cleaned
    fi
    echo "You can find your cleaned fastq files at $out/cleaned" 
fi

# Create directory to store trimmed fastq files if pretrim is not specified
if [ ! "$pretrimmed_flag" = 'T' ]; then
    if [ ! -d "$out/trimmed" ] 
    then
        mkdir $out/trimmed
    fi
    echo "You can find your trimmed fastq files at $out/trimmed" 
fi

# Create directory to store fastqc results
if [ ! -d "$out/fastqc_results" ] 
then 
    mkdir $out/fastqc_results
fi
echo "You can find all the fastqc test results at $out/fastqc_results"

# Main loop
while IFS= read -r line
do
    echo "begin downloading $line..."
    $dir/bin/fasterq-dump "$line" -O $out --verbose

    # Case where files are not paired-end, needs to be cleaned first
    if [ ! "$paired_end_flag" = 'PE' ]; then

        echo "combining $line into a single fastq file..."
        cat $out/$line*.fastq > $out/cleaned/$line.fastq

        echo "performing fastqc on combined file..."
        fastqc $out/cleaned/$line.fastq -o $out/fastqc_results

        if [ ! "$pretrimmed_flag" = 'T' ]; then
        echo "trimming..."
        ~/TrimGalore-0.6.10/trim_galore -o $out/trimmed $out/cleaned/$line.fastq

        echo "trimming complete, performing fastqc..."
        fastqc $out/trimmed/$line"_trimmed.fq" -o $out/fastqc_results
        fi

    # Case where files are paired end, no cleaning needed
    else 
        echo "performing fastqc on fastq file..."
        fastqc $out/$line.fastq -o $out/fastqc_results 

        if [ ! "$pretrimmed_flag" = 'T' ]; then
        echo "trimming..."
        ~/TrimGalore-0.6.10/trim_galore -o $out/trimmed $out/$line.fastq

        echo "trimming complete, performing fastqc..."
        fastqc $out/trimmed/$line"_trimmed.fq" -o $out/fastqc_results
        fi

    fi
    echo "$line complete"
    echo ""
done < "$txt"
