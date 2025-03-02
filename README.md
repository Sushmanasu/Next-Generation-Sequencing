# Next-Generation-Sequencing

NGS_WES_Pipeline

Genomic Data Analysis Pipeline

Introduction

This repository provides a comprehensive pipeline for genomic data analysis, specifically focusing on whole-genome
sequencing and variant calling using the Genome Analysis Toolkit (GATK). It includes steps for Quality control, read trimming, genome assembly or reference mapping, variant calling, and annotation. This pipeline utilizes Docker for GATK ensure consistent and reproducible results.

The provided commands and scripts guide users through the entire workfiow from raw FASTQ filies to annotated variant calls, using widiely adopted bioinformatics tools and best practices.

Requirements


To run this pipeline, ensure the following are installed:


* Docker: For running the GATK container and other tools in isolated environments.
* GATK: The Genome Analysis Toolkit for variant discovery in high-throughput sequencing data.
* BWA: For aligning short DNA sequences to a reference genome.
* SAMtools: For manipulating SAM/BAM files
* FASTQC: For quality control of raw sequencing data
* FASTP: For quality trimming and adapter removal from FASTQ files
* ANNOVAR: For functional annotation of genetic variants.

Pipeline

Step 1: Download the SARS-CoV-2 Sample

•	Explanation: The wget command is used to download raw sequencing data in FASTQ format from the European Nucleotide Archive (ENA).
________________________________________
Step 2: Install FastQC

•	Explanation: FastQC is a tool for checking the quality of sequencing reads. 
________________________________________
Step 3: Alternative Method – Download Using fastq-dump

•	Explanation: The sra-toolkit provides an alternative way to download sequencing data directly from the NCBI Sequence Read Archive (SRA).
________________________________________
Step 4: Download the Reference Genome

•	Explanation: Downloads the reference genome from NCBI in FASTA format.
________________________________________
Step 5: Run Quality Check on FASTQ Files Using FastQC 

•	Explanation: Runs FastQC on both paired-end sequencing files to assess sequence quality, GC content, adapter contamination, and other quality metrics.
________________________________________
Step 6: Understand the FASTQ and FASTQC Reports

•	FASTQ format: Contains raw sequencing reads and their quality scores.

•	FASTQC report: Provides visual and statistical summaries of sequence quality, including per-base sequence quality and overrepresented sequences.
________________________________________
Step 7: Mapping (Create a Directory for Mapping)

•	Explanation: Creates a directory named Mapping where alignment results will be stored.
________________________________________
Step 8: Index the Reference Genome and Perform Alignment

•	Explanation:

o	bwa index: Creates an index of the reference genome to allow efficient alignment.
o	bwa mem: Aligns the sequencing reads to the reference genome and outputs a Sequence Alignment Map (SAM) file.
________________________________________
Step 9: Convert, Sort, and Index the Alignment

samtools view -@ 20 -S -b Mapping/ERR5743893.sam > Mapping/ERR5743893.bam

•	Explanation: Converts the SAM file into a BAM (binary alignment map) format using samtools view. The -@ 20 flag uses 20 CPU threads for faster processing.

samtools sort -@ 32 -o Mapping/ERR5743893.sorted.bam Mapping/ERR5743893.bam

•	Explanation: Sorts the BAM file for efficient processing in downstream analyses.

samtools index Mapping/ERR5743893.sorted.bam

•	Explanation: Indexes the sorted BAM file for fast access in visualization and variant calling.
________________________________________
Step 10: Visualization in IGV

•	Explanation: The sorted BAM file can be loaded into Integrative Genomics Viewer (IGV) to visually inspect the read alignments and detect potential issues such as coverage gaps or misalignments.
________________________________________
Step 11: Index the Reference Genome

•	Explanation: Creates an index for the reference genome FASTA file to speed up querying specific genomic regions.
_______________________________________
Step 12: Variant Analysis

freebayes 

•	Explanation: Runs FreeBayes, a Bayesian variant caller, to identify SNPs and small indels in the sequencing data compared to the reference genome. The output is in Variant Call Format (VCF).

•	Explanation:

	bgzip: Compresses the VCF file.
  tabix: Indexes the compressed VCF file for rapid querying.
________________________________________
Final Outcome

By following these steps, you will:

•	Download sequencing reads and reference genome

•	Perform quality control

•	Align reads to the reference genome

•	Sort and index the alignments

•	Visualize the mapping results in IGV

•	Call genetic variants from the sequencing data

