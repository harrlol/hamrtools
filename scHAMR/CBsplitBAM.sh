#!/bin/bash
set -eu

# This script takes in the output of 10X Genomic's cell ranger count command and splits it into clusters, as determined by the upstream
# R analyses. This prepares the bam for more fruitful HAMR analysis.

if [ "$#" -lt 4 ]; then
echo "Missing arguments!"
echo "USAGE: CBsplitBAM.sh <possorted_genome_bam.bam> <CBs_Clusters_dataframe.csv> <CB dir> <cellranger path>"
echo "EXAMPLE:"
exit 1
fi

bam=$1
csv=$2
csvfn=$(basename "$csv")
out=$3
cr=$4

# Create necessary directories
if [ ! -d $out/CB ] 
then
    mkdir $out/CB
    echo "created path: $out/CB"
fi

if [ ! -d $out/CB/raw ] 
then
    mkdir $out/CB/raw
    echo "created path: $out/CB/raw"
fi

if [ ! -d $out/CB/csv ] 
then
    mkdir $out/CB/csv
    echo "created path: $out/CB/csv"
fi

if [ ! -d $out/CB/txt ] 
then
    mkdir $out/CB/txt
    echo "created path: $out/CB/txt"
fi

# Create CB divided bam file folder
if [ ! -d $out/CB/split_bam_out ] 
then
    mkdir $out/CB/split_bam_out
    echo "created path: $out/CB/split_bam_out"
fi

# Create temp file folder
if [ ! -d $out/CB/split_temp ] 
then
    mkdir $out/CB/split_temp
    echo "created path: $out/CB/split_temp"
fi

# Extract info from csv
cp $csv $out/CB/raw
sed -i 1d $out/CB/raw/$csvfn  #remove header 
cd $out/CB/raw
awk '{gsub(/"/, "", $2); print > ("c" $2 ".csv")}' FS=, $out/CB/raw/$csvfn
echo "extracting info from csv file..."
rm $out/CB/raw/$csvfn

# Populate ~/csv/ with CB from .csv in ~/raw/ and split each cluster into a single csv, while also counting 
cd $out/CB/raw
nClu=1
for i in $(ls *.csv)
do 
  cut -d\  -f1  $i > $out/CB/csv/${i}
  ((nClu++))
done
nClu=$((nClu - 2))
cd

echo "formatting each csv cluster to txt..."
# Populate ~/txt/ by adding CB:Z to every line of each .csv and remove double quotes
for i in $(seq 0 $nClu)
do 
  sed -e 's/^/CB:Z:/' -e 's/"//g' $out/CB/csv/c${i}.csv > $out/CB/txt/c${i}.txt
done

# The 10x bam output (use unfiltered)
export BAM_FILE=$bam

cd $out/CB/txt
for i in $(ls *.txt)
do
    filename_ext=$(basename "$i")
    filename=${filename_ext%.*}
    echo "splitting bam files for cluster: $filename..."
    
    # Save the header lines
    samtools view -H $BAM_FILE > $out/CB/split_temp/${filename}_SAM_header

    # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
    samtools view $BAM_FILE | LC_ALL=C grep -F -f c13.txt > $out/CB/split_temp/${filename}_filtered_SAM_body

    # Combine header and body
    cat $out/CB/split_temp/${filename}_SAM_header $out/CB/split_temp/${filename}_filtered_SAM_body > $out/CB/split_temp/${filename}_filtered.sam

    # Convert filtered.sam to BAM format
    samtools view -b $out/CB/split_temp/${filename}_filtered.sam > $out/CB/split_bam_out/${filename}_filtered.bam
done
