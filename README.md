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
