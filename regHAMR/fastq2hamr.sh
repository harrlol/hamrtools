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
echo "USAGE: fastq2hamr.sh <sample.fastq> <anotation.gtf/gff3> <bowtie dir> 
<genome.fasta> <sorting_script.pl> <hamr_model.Rdata> <out dir> <mismatch_num> <filename dict>"
echo "EXAMPLE:"
exit 1
fi

smp=$1
ant=$2
bt=$3
gno=$4
sort=$5
mdl=$6
#root out, true out and hamr out assigned later
rout=$7
len=$8
dic=$9
# Assume single end
det=1

smpext=$(basename "$smp")
smpdir=$(dirname "$smp")
smpkey="${smpext%.*}"
smpname=""
original_ext="${smpext##*.}"

if [[ $smpkey == *_1* ]]; then
  smpkey="${smpkey%_1*}"
  smp1="$smpdir/${smpkey}_1.$original_ext"
  smp2="$smpdir/${smpkey}_2.$original_ext"
  # Paired end recognized
  det=0
  echo "$smpext is a part of a paired-end sequencing file"
elif [[ $smpkey == *_2* ]]; then
  # If _2 is in the filename, this file was processed along with its corresponding _1 so we skip
  echo "$smpext has already been processed with its _1 counter part. Skipped."
  echo ""
  exit 1
else 
  echo "$smpext is a single-end sequencing file"
  echo ""
fi

# Real the CSV file into a DataFrame
mapfile -t names < <(awk -F, '{ print $1 }' "$dic")
mapfile -t smpf < <(awk -F, '{ print $2 }' "$dic")

# Create a dictionary from the DataFrame
declare -A dictionary
for ((i=0; i<${#names[@]}; i++)); do
    dictionary[${names[i]}]=${smpf[i]}
done


if [[ $smpkey == *_trimmed* ]]; then
  smpkey="${smpkey%_trimmed*}"
fi

# Retrieve the translated value
if [[ ${dictionary[$smpkey]+_} ]]; then
    smpname="${dictionary[$smpkey]}"
    smpname="${smpname//$'\r'}"
    echo "Sample Group: $smpname"
else
    echo "No filename found for input: $smpkey"
    exit 1
fi

# Reassign / declare pipeline file directory
if [ ! -d "$rout/pipeline/$smpkey""_temp" ] 
then
    mkdir "$rout/pipeline/$smpkey""_temp"
    echo "created path: $rout/pipeline/$smpkey""_temp"
fi

out=$rout/pipeline/$smpkey"_temp"
echo "You can find all the intermediate files for $smpkey at $out" 


# Reassign hamr output directory
if [ ! -d "$rout/hamr_out" ] 
then
    mkdir $rout/hamr_out
    echo "created path: $rout/hamr_out"
fi

hamrout=$rout/hamr_out
echo "You can find the HAMR output file for $smpkey at $hamrout/$smpname.mod.txt" 


echo ""
echo "sample path = $1"
echo "sample SRR = $smpkey"
echo "sample name = $smpname"
echo "annotation = $2"
echo "bowtie directory = $3"
echo "genome = $4"
echo "scorting script = $5"
echo "hamr model = $6"
echo "mismatch number = $8"
echo ""


red=8
if (($len > 8)); then
  red=$((len +1))
fi

#maps the trimmed reads to provided annotated genome, can take ~1.5hr
#fine tuning of the arguments might be needed
if [ "$det" -eq 1 ]; then
  echo "Performing tophat2 with a single-end file."
  tophat2 --library-type fr-firststrand --read-mismatches $len --read-edit-dist $red  --max-multihits 10 \
    --b2-very-sensitive \
    --transcriptome-max-hits 10 \
    --no-coverage-search \
    --output-dir $out \
    -G $ant \
    -p 4 \
    $bt \
    $smp

else
  echo "Performing tophat2 with a paired-end file."
  tophat2 --library-type fr-firststrand --read-mismatches $len --read-edit-dist $red --max-multihits 10 \
    --b2-very-sensitive \
    --transcriptome-max-hits 10 \
    --no-coverage-search \
    --output-dir $out \
    -G $ant \
    -p 4 \
    $bt \
    $smp1 $smp2
fi


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
    -fe $out/unique_RG_ordered_splitN.resort.bam $gno $mdl $out $smpname 30 50 0.01 H4 .01 0 .05

if [ ! -e "$out/${smpname}.mods.txt" ]
then 
    cd $hamrout
    printf "${smpname} \n" >> zero_mod.txt
    cd
else
# HAMR needs separate folders to store temp for each sample, so we move at the end
    cp $out/${smpname}.mods.txt $hamrout
fi
