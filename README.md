# Development of metagenomics tools to deciphering the virome and genomic diversity in natural isolates

### By: Marta Ibañez Lligoña

### Supervised by: Dr Josep Quer

Universitat Autònoma de Barcelona (UAB) and Vall d'Hebron Institut de Recerca (VHIR). Grup Malalties Hepàtiques-Hepatitis viral.

### More information: missing information and results in this repository will be kept private until the paper is published. 
## Repository content

### Readme: 
Information about repository and folders content.

## Scripts ##
### Bioinformatics pipeline ###
This folder contains all scripts used for the analysis of fastQ files from metagenomics samples (Culex sp. & G. affinis).

| Script file | Description | 
| ----------------------------- | ----------------------- | 
trimmomatic.run | Script used to remove low quality, short reads and overrepresented sequences from fastQ files. | 
bowtie2-build.CP.run | Script used to build C. pipens bowtie2 database. |
bowtie2-build.CQ.run | Script used to build C. quinquefasciatus bowtie2 database. | 
bowtie2-build.GA.run | Script used to build G. affinis bowtie2 database. |
bowtie2-buildHS.run | Script used to build human bowtie2 database. |
bowtie2-mapCP.run | Script used to map reads to C. pipiens bowtie2 database . |
bowtie2-mapCQ.run | Script used to map reads to C. quinquefasciatus bowtie2 database. |
bowtie2-mapGA.run | Script used to map reads to G. affinis bowtie2 database. |
bowtie2-mapHS.run | Script used to map reads to human bowtie2 database. |
fastq_screen.conf | Configuration file to run fastQ-screen in mosquito samples. |
fastq_screen.run | Script used to run fastQ-screen and identify host genomes in mosquito and fish samples. |
fastq_screenGA.conf |Configuration file to run fastQ-screen in fish samples. |
samtools_bam2fastq_sort.run | Script used to sort bam file. |
samtools_bam2fastq.run | Script used to transform bam file to separated forward and reverse fastQ files. |
samtools_mapping.run | Script used to transform sam file to bam. |
samtools_unmapped_reads.run | Script used to remove reads mapping to host genome. |

## Interpretation of results ##
This folder contains all scripts used for analysis of results and making of graphs.

| Script file | Description | 
| ----------------------------- | ----------------------- | 
kaiju_analysis.R | Script used to analyse results from kaiju of main metagenomic analysis from C. pipiens and G. affinis. | 
kaiju_graphs.R | Script used to analyse results from kaiju of main metagenomic analysis from C. pipiens and G. affinis. | 
QA_control.R | Script used to analyse results from from the sensitivity analysis. | 
pipeline_graphs.R | Scripts used to make graphs to show effectiveness of the bioinformatics pipeline. |