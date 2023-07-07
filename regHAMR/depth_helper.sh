#!/bin/bash
set -u

# Script to use sratool (newest version) that is otherwise unavailable in bash 

# First, install the relative sratoolkit manually from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

# Next, create an SRR accession list (in .txt format) using the SRA run selector tool.

# Locate both your sratoolkit and your .txt


if [ "$#" -lt 5 ]; then
echo "Missing arguments!"
echo "EXAMPLE:"
exit 1
fi

out=$1
t=$2
d=$3
n=$4
f=$5

for ff in $out/pipeline/depth/*.bam
    do
        if echo "$ff" | grep -q "$n"
        then
            tt=$(basename $ff)
            nn=${tt%.*}
            echo "[$n] extracting depth information from $nn"
            for i in $(seq 1 $(wc -l < $f))
            do
                chr=$(sed "${i}q;d" $f | sed 's/\t/\n/g' | sed '1q;d')
                pos=$(sed "${i}q;d" $f | sed 's/\t/\n/g' | sed '2q;d')
                dph=$(samtools coverage \
                    -r $chr:$pos-$pos \
                    $ff \
                    | awk 'NR==2' | awk -F'\t' '{print $7}')
                awk -v "i=$i" 'NR==i {print $0"\t"var; next} 1' var="$dph" $f > $d/${nn}_new.bed && mv $d/${nn}_new.bed $f &
            done
            echo "[$n] finished $nn"
        fi
    done