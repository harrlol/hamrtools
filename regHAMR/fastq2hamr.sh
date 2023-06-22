#!/bin/bash
set -u

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

if [ "$#" -lt 11 ]; then
echo "Missing arguments!"
echo "USAGE: fastq2hamr.sh <sample.fastq> <anotation.gtf/gff3> 
<genome.fasta> <sorting_script.pl> <hamr_model.Rdata> <out dir> <mismatch_num> <filename dict> <thread num> <genome length> <repo>"
echo "EXAMPLE:"
exit 1
fi

smp=$1
ant=$2
gno=$3
sort=$4
mdl=$5
#root out, true out and hamr out assigned later
rout=$6
mis=$7
ohang=$((mis-1))
dic=$8
cores=$9
gnolen=${10}
repo=${11}
# Assume single end
det=1

#######################
# Running checks to ensure program can run normally
gatk=""
hamr=""

# Grab command from repo where needed
shopt -s nocaseglob

for folder in "$repo"/*; do
    # Check if the folder name contains the string "sratoolkit"
    if [[ "$folder" == *gatk* ]]; then
        gatk="$folder"
    # Check if the folder name contains the string "hamr_py3" (case unsensitive)
    elif [[ "${folder,,}" == *hamr_py3* ]]; then
        hamr="$folder"
    fi
done

shopt -u nocaseglob

if [ ! -n "$gatk" ]; then
    echo "GATK not installed, please check."
    exit 1
fi

if [ ! -n "$hamr" ]; then
    echo "HAMR not installed, please check."
    exit 1
fi


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

# Read the CSV file into a DataFrame
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
    echo "[$smpkey] Sample group name found: $smpname"
else
    echo "[$smpkey] Could not locate sample group name, exiting..."
    exit 1
fi

# Reassign / declare pipeline file directory
if [ ! -d "$rout/pipeline/$smpkey""_temp" ] 
then
    mkdir "$rout/pipeline/$smpkey""_temp"
    echo "[$smpkey] created path: $rout/pipeline/$smpkey""_temp"
fi

out=$rout/pipeline/$smpkey"_temp"
echo "[$smpkey] You can find all the intermediate files for $smpkey at $out" 


# Reassign hamr output directory
if [ ! -d "$rout/hamr_out" ] 
then
    mkdir $rout/hamr_out
    echo "created path: $rout/hamr_out"
fi

hamrout=$rout/hamr_out
echo "[$smpkey] You can find the HAMR output file for $smpkey at $hamrout/$smpname.mod.txt" 


echo "[$smpkey] Begin HAMR pipeline"
cd "$rout/pipeline/$smpkey""_temp"
#maps the trimmed reads to provided annotated genome, can take ~1.5hr
#fine tuning of the arguments might be needed
if [ "$det" -eq 1 ]; then
  echo "[$smpkey] Performing STAR with a single-end file."
  STAR \
    --runThreadN $cores \
    --genomeDir $rout/ref \
    --readFilesIn $smp \
    --sjdbOverhang $ohang \
    --sjdbGTFfile $ant \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax $mis \
    --outSAMtype BAM SortedByCoordinate
else
  echo "[$smpkey] Performing STAR with a paired-end file."
  STAR \
    --runThreadN $cores \
    --genomeDir $rout/ref \
    --readFilesIn $smp1 $smp2 \
    --sjdbOverhang $ohang \
    --sjdbGTFfile $ant \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax $mis \
    --outSAMtype BAM SortedByCoordinate
fi
cd

wait

#sorts the accepted hits
echo "[$smpkey] sorting..."
samtools sort \
    -n $out/Aligned.sortedByCoord.out.bam \
    -o $out/sort_accepted.bam
echo "[$smpkey] finished sorting"
echo ""

wait

#filter the accepted hits by uniqueness
echo "[$smpkey] filter unique..."
samtools view \
    -h $out/sort_accepted.bam \
    | perl $sort 1 \
    | samtools view -bS - \
    | samtools sort \
    -o $out/unique.bam
echo "[$smpkey] finished filtering"
echo ""

wait

#adds read groups using picard, note the RG arguments are disregarded here
echo "[$smpkey] adding/replacing read groups..."
gatk AddOrReplaceReadGroups \
    I=$out/unique.bam \
    O=$out/unique_RG.bam \
    RGID=1 \
    RGLB=xxx \
    RGPL=illumina_100se \
    RGPU=HWI-ST1395:97:d29b4acxx:8 \
    RGSM=sample
echo "[$smpkey] finished adding/replacing read groups"
echo ""

wait

if [ ! -e "$rout/ref/picard_ref.dict" ] 
then 
    gatk CreateSequenceDictionary \
    R=$gno \
    O=$rout/ref/picard_ref.dict
    echo "Picard genome dictionary created"
fi

wait

#reorder the reads using picard
echo "[$smpkey] reordering..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=$out/tmp" ReorderSam \
    I=$out/unique_RG.bam \
    O=$out/unique_RG_ordered.bam \
    R=$gno \
    CREATE_INDEX=TRUE \
    SEQUENCE_DICTIONARY=$rout/ref/picard_ref.dict \
    TMP_DIR=$out/tmp
echo "[$smpkey] finished reordering"
echo ""

wait

#splitting and cigarring the reads, using genome analysis tool kit
#note can alter arguments to allow cigar reads 
echo "[$smpkey] getting split and cigar reads..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=$out/tmp" SplitNCigarReads \
    -R $gno \
    -I $out/unique_RG_ordered.bam \
    -O $out/unique_RG_ordered_splitN.bam \
    # -U ALLOW_N_CIGAR_READS
echo "[$smpkey] finished splitting N cigarring"
echo ""

wait

#final resorting using picard
echo "[$smpkey] resorting..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=$out/tmp" SortSam \
    I=$out/unique_RG_ordered_splitN.bam \
    O=$out/unique_RG_ordered_splitN.resort.bam \
    SORT_ORDER=coordinate
echo "[$smpkey] finished resorting"
echo ""

wait

#hamr step, can take ~1hr
echo "[$smpkey] hamr..."
python $hamr/hamr.py \
    -fe $out/unique_RG_ordered_splitN.resort.bam $gno $mdl $out $smpname 30 50 0.01 H4 1 .05 .05

wait

if [ ! -e "$out/${smpname}.mods.txt" ]
then 
    cd $hamrout
    printf "${smpname} \n" >> zero_mod.txt
    cd
else
# HAMR needs separate folders to store temp for each sample, so we move at the end
    cp $out/${smpname}.mods.txt $hamrout
fi

# Move the unique_RG_ordered.bam and unique_RG_ordered.bai to a folder for read depth analysis
cp $out/unique_RG_ordered.bam $rout/pipeline/depth/$smpname.bam
cp $out/unique_RG_ordered.bai $rout/pipeline/depth/$smpname.bai
