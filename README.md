# 5505_Group_Project: RNA Sequencing and Differential Abundance Analysis Pipeline
### Alina Elahie, Aleksandra Maciruta, Kaleem Uddin Mohammed, Rachel Yorke
## Overview
This project analyzes RNA-Seq data from **GSE270472** to study transcriptional deregulation in Huntington's Disease (HD). The pipeline processes sequencing data, identifies differentially expressed genes (DEGs), and performs pathway analysis.

## Workflow
1. **Quality Control**: FastQC is used to check read quality.
2. **Trimming**: BBDUK 2 removes low-quality bases and adapter sequences.
3. **Alignment**: Bowtie2 aligns reads to the GRCh38 reference genome.
4. **Quantification**: RSEM estimates transcript abundances.
5. **Differential Expression Analysis**: DESeq2 is used to find DEGs.
6. **Pathway Analysis**: GSEA and gProfiler2 identify enriched pathways and GO terms.

## Reference Genome & Annotations
- **Genome:** GRCh38 (ENSEMBL v102)
- **Annotation:** ENSEMBL GTF file

## Software & Tools
- **Trimming:** BBDUK 2
- **Alignment:** Bowtie2
- **Quantification:** RSEM
- **Differential Analysis:** DESeq2
- **Pathway Analysis:** GSEA, gProfiler2

## Usage
Run the **Nextflow pipeline** using:
```bash
nextflow run main.nf --reads 'path/to/reads/*_R{1,2}.fastq.gz' --genome 'path/to/genome.fa' --outdir 'path/to/output/'
