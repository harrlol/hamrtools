# hamrbox
a box of tools that might be needed before or after HAMR

## fastq2hamr
a shell wrapper program that takes in an adaptor-trimmed fastq file all the way to applying HAMR, yielding intermediates and hamr prediction table
this is representative of the operational pipeline used to analyze data generated from wt-mta GMUCT and RNA-seq. The pipeline includes the preparing/cleaning steps recommended by Vandivier et al. https://link.springer.com/protocol/10.1007/978-1-4939-8808-2_4

USAGE: 
```
fastq2hamr.sh <sample.fastq> <anotation.gtf/gff3> <bowtie dir> <genome.fasta> <sorting_script.pl> <hamr_model.Rdata> <out dir> <mismatch_num>
```
```sample.fastq```: a fastq file that is adaptor-processed using trim-galore or cutadapt. Use fastqc to monitor the quality of trimmed and pre-trimmed fastq files in case of anomalies.

```anotation.gtf/gff3```: these are specific to the model organism you are working with and can be obtained through various databanks, a simple googling should work.

```bowtie dir```: this requires you to build a bowtie dictionary before-hand, via:
```
bowtie2-build <genome.fa> <output directory>
```
In addition, a dictionary is needed along some steps, create via:
```
java YOURDIR/picard.jar CreateSequenceDictionary R=genome.fa O=genome.dict
```

```genome.fasta```: also specific to your model organism, can be obtained via googling.

```sorting_script.pl```: a script developed for internal use.

```hamr_model.Rdata```: can be found under HAMR/models/euk_trna_mods.Rdata.

```out dir```: the directory where your intermediate files and your output table will be located.

```mismatch_num```: gather the sequencing length from the fastqc step, use 0.06 to multiply the sequencing length to obtain this number

## consensusOverlap
Script to overlap consensus mods with various libraries (UTR, CDS, gene, etc.) after post HAMR processing in R, the libraries are obtained in various banks like TAIR.

USAGE: 
```
consensusOverlap.sh <consensus.bed> <cds.bed> <utr.bed> <gene.bed> <mrna.bed> <out dir>
```
```consensus.bed```: a bed file obtained from processing hamr predicted mod tables and keeping union/intersection of mods between biological replicates. 

```cds.bed | utr.bed | gene.bed | mrna.bed```: annotated coding region for a given organism, obtainable via gtf/gff3 files. Please see https://github.com/dtrain16/NGS-scripts/blob/master/TAIR10_annotation.sh for example.

```out dir```: where you want your overlap outputs to be, note that the bed file will be individually overlapped with all 4 library bed annotations, so 4 outputs will be generated. 

## fasterq-dumpAdaptor
We had difficulty installing sratoolkit and using the newest functions, so I created this script to manually use it but still high throughput. A quick tutorial:

First, install the relative sratoolkit manually from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

Next, create an SRR accession list (in .txt format) using the SRA run selector tool. Locate both your sratoolkit and your .txt.

USAGE: 
```
fasterq-dumpAdaptor.sh <sratoolkit dir> <accession.txt> <out dir>
```
```sratoolkit dir```: where your sratoolkit is located.
```accession.txt```: where your accession txt list is located.
```out dir```: where you want your final fastq files to be.

