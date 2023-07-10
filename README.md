# hamrbox
This is a toolbox tailored for the High Throughput Annotation of Modified Ribonucleotides, abbreviated HAMR, developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684).
The main function of hamrbox is to make the original method more accessible by automating the tedious pre-processing steps, allowing 
users to analyze entire experiments at a time. More importantly, this repository adapts the HAMR pipeline to single-cell sequencing data, 
while integrating Seurat as a part of the clustering analysis. 

This tool aims to take in as simple as an SRR accession code, and returns HAMR analysis in highly visualized, intuitive plots or tables that aid in downstream analysis or biological discoveries. 

hamrbox is split into the below 2 suites: 
### regHAMR
regular HAMR, crystalizes all the steps in Ryvkin et al.'s paper and extends the analysis through a dozen algorhithms in R
that transforms the HAMR output in more accessible ways.

### scHAMR
single-cell HAMR, includes crucial preprocessing steps neccesitated for single cell data types. This method uses cell ranger
to align the reads, Seurat to cluster the reads by cell type, then HAMR the reads by clusters. 

Please see the readme inside each suite for a more detailed instruction on how to use the respective programs. 
