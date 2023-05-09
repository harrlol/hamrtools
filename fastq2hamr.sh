#!/bin/bash
set -eu

# Script to process sequencing data in fastq format to apply HAMR. See 
#HAMR: High-Throughput Annotation of Modified Ribonucleotides, Lee E. Vandivier et al. for reference. 

### CONDA environment is installed
# create --name py2
# conda install -n py2 -c bioconda gatk
# conda install -n py2 -c bioconda bowtie2
# conda install -n py2 -c bioconda tophat2
# conda install -n py2 -c bioconda bedtools
# conda install -n py2 -c bioconda conda-forge/python
# conda install -n py2 -c conda-forge python2.7
# conda install -n py2 r-class
# conda install -n py2 -c bioconda subread
# conda install -n py2 -c bioconda bioconvert
# conda install -n py2 -c r r-tidyverse
# conda activate py2

#note bowtie index needs to be prebuilt and directory available prior to running this script
#input fastq file should be pre-trimmed of adaptors

if [ "$#" -lt 8 ]; then
echo "Missing arguments!"
echo "USAGE: fastq2hamr.sh <sample.fastq> <anotation.gtf/gff3> <bowtie dir> <genome.fasta> <sorting_script.pl> <hamr_model.Rdata> <out dir> <read_length>"
echo "EXAMPLE:"
exit 1
fi

smp=$1
ant=$2
bt=$3
gno=$4
sort=$5
mdl=$6
out=$7
len=$8

echo ""
echo "sample = $1"
echo "annotation = $2"
echo "bowtie directory = $3"
echo "genome = $4"
echo "scorting script = $5"
echo "hamr model = $6"
echo "out = $7"
echo "read length = $8"
echo ""

#maps the trimmed reads to provided annotated genome, can take ~1.5hr
#fine tuning of the arguments might be needed
echo "mapping reads to genome using tophat2..."
tophat2 --library-type fr-firststrand --read-mismatches $len --read-edit-dist 8 --max-multihits 10 \
    --b2-very-sensitive \
    --transcriptome-max-hits 10 \
    --no-coverage-search \
    --output-dir $out \
    -G $ant \
    -p 4 \
    $bt \
    $smp

#sorts the accepted hits
echo "sorting..."
samtools sort \
    -n $out/accepted_hits.bam \
    -o $out/sort_accepted.bam
echo "finished sorting"
echo ""

#filter the accepted hits by uniqueness
echo "filter unique..."
samtools view \
    -h $out/sort_accepted.bam \
    | perl $sort 1 \
    | samtools view -bS - \
    | samtools sort \
    -o $out/unique.bam
echo "finished filtering"
echo ""

#adds read groups using picard, note the RG arguments are disregarded here
echo "adding/replacing read groups..."
java -jar /usr/bin/picard/picard.jar AddOrReplaceReadGroups \
    I=$out/unique.bam \
    O=$out/unique_RG.bam \
    RGID=1 \
    RGLB=xxx \
    RGPL=illumina_100se \
    RGPU=HWI-ST1395:97:d29b4acxx:8 \
    RGSM=sample
echo "finished adding/replacing read groups"
echo ""

#reorder the reads using picard
echo "reordering..."
java -Xmx2g -Djava.io.tmpdir=$out/tmp \
    -jar /usr/bin/picard/picard.jar ReorderSam \
    I=$out/unique_RG.bam \
    O=$out/unique_RG_ordered.bam \
    R=$gno \
    CREATE_INDEX=TRUE \
    TMP_DIR=$out/tmp
echo "finished reordering"
echo ""

#splitting and cigarring the reads, using genome analysis tool kit
#note can alter arguments to allow cigar reads 
echo "getting split and cigar reads..."
/Data04/harrli02/repo/gatk-4.3.0.0/gatk SplitNCigarReads \
    -R $gno \
    -I $out/unique_RG_ordered.bam \
    -O $out/unique_RG_ordered_splitN.bam
    #  -U ALLOW＃_N_CIGAR_READS
echo "finished splitting N cigarring"
echo ""

#final resorting using picard
echo "resorting..."
java -Xmx2g -Djava.io.tmpdir=$out/tmp  \
    -jar /usr/bin/picard/picard.jar SortSam \
    I=$out/unique_RG_ordered_splitN.bam \
    O=$out/unique_RG_ordered_splitN.resort.bam \
    SORT_ORDER=coordinate
echo "finished resorting"
echo ""

#hamr step, can take ~1hr
echo "hamr..."
python2.7 /Data04/harrli02/repo/HAMR/hamr.py \
    -fe $out/unique_RG_ordered_splitN.resort.bam $gno $mdl $out hamr_out 30 1 0.01 H4 1 .05 .05

fi