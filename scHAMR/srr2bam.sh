#!/bin/bash
set -eu

# This script takes in an SRR accession code, and sets up the necessary paths/directories in the folder
# assigned, and populate those folders with necessary intermediate files for seurat analysis.

if [ "$#" -lt 6 ]; then
echo "Missing arguments!"
echo "USAGE: srr2bam.sh <SRR Accession> <sra toolkit dir> <cr dir> <genome.fa> <annotation.gtf> <outdir>"
echo "EXAMPLE:"
exit 1
fi

srr=$1
srakit=$2
cr=$3
gno=$4
ant=$5
out=$6

if [ ! -d $out ] 
then
    mkdir $out
    echo "created path: $out"
fi

if [ ! -d $out/datasets ] 
then
    mkdir $out/datasets
    echo "created path: $out/datasets"
fi

cd $out/datasets
# Locate the dataset in question, and fasterq dump
$srakit/bin/fasterq-dump \
    $srr \
    --split-files \
    --include-technical \
    --verbose
cd

if [ ! -d $out/pipeline ] 
then
    mkdir $out/pipeline
    echo "created path: $out/pipeline"
fi

if [ ! -d $out/pipeline/ref ] 
then
    mkdir $out/pipeline/ref
    echo "created path: $out/pipeline/ref"
fi

cp $ant $out/pipeline/ref
cp $gno $out/pipeline/ref
ant=$out/pipeline/ref/$(basename "$ant")
gno=$out/pipeline/ref/$(basename "$gno")
antf="${ant%.*}_filtered.${ant##*.}"

# Filter gtf for protein coding genes
$cr/bin/cellranger mkgtf \
    $ant \
    $antf \
    --attribute=gene_biotype:protein_coding
    
# Make ref files for count
cd $out/pipeline/ref
$cr/bin/cellranger mkref \
    --genome=reannotated \
    --fasta=$gno \
    --genes=$antf

cd $out/datasets
for i in $(ls *.fastq)
do
    t=${i%.*}
    v=${t##*_}
    s=${i%%_*}
    if [ "$v" -lt 3 ]; then
        mv $i ${s}_S1_L001_R${v}_001.fastq
    else
        mv $i ${s}_S1_L001_I$((v-2))_001.fastq
    fi
done

# Align reads to genome using cell ranger (slow, 2hr)
cd $out/pipeline
$cr/bin/cellranger count --id=cell_ranger_out \
    --transcriptome=$out/pipeline/ref/reannotated \
    --fastqs=$out/datasets \
    --sample=$srr \
    --localcores=8 \
    --localmem=64
