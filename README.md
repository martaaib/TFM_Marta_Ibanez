# Identification of infectious agent genomes and implementation of a sensitivity study in Culex sp. and G. affinis 

### By: Marta Ibañez Lligoña

### Supervised by: Dr Josep Quer

Universitat Autònoma de Barcelona (UAB) and Vall d'Hebron Institut de Recerca (VHIR). Grup Malalties Hepàtiques.

### More information: missing information and results in this repository will be kept private until the paper is published. 
## Repository content

### Readme: 
Information about repository and folders content.

## Scripts ##
### Bioinformatics pipeline ###
| Script file | Description | 
| ----------------------------- | ----------------------- | 
trimmomatic.run | Script used to remove low quality, short reads and overrepresented sequences from fastQ files. | 
bowtie2-build.CP.run | Script used to build C. pipens bowtie2 database. |
bowtie2-build.CQ.run | Script used to build C. quinquefasciatus bowtie2 database. | 
bowtie2-build.GA.run | Script used to build G. affinis bowtie2 database. |
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
