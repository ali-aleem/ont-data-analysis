# 🧬 Human Small Variants Pipeline (WSL)
This repository provides a reproducible workflow for **Single Nucleotide Variants (SNVs)** and **Insertions/Deletions (InDels)** from **Oxford Nanopore Technologies (ONT)** raw reads.  
The pipeline covers **preprocessing**, **quality control**, **alignment**, **variant calling**, and **annotation** steps, suitable for both beginners and advanced users.

## 📑 Table of Contents
1. [Overview & Workflow](#1-overview--workflow)
2. [Pipeline Tools](#2-pipeline-tools)
3. [Installation Guide](#3-installation-guide)
   - [A. WSL Setup](#a-wsl-setup)
   - [B. Miniconda Setup](#b-miniconda-setup)
   - [C. Tools Installation](#c-tools-installation)
   - [D. Reference Genome and Databases](#d-reference-genome-and-databases)
4. [Pipeline Commands](#4-pipeline-commands)
5. [References](#5-references)

## 1. Overview & Workflow
This pipeline performs **accurate and efficient detection of SNVs and InDels** from ONT sequencing data.  
It runs seamlessly on **Windows (via WSL)** and standard Linux environments.
### **Input Files**
- Raw ONT FASTQ reads (`*.fastq.gz`)
- Reference genome (`.fa`)
### **Workflow Summary**

| Step | Description | Tools |
|------|--------------|-------|
| **1. Data Preprocessing** | Merge barcoded reads | `cat` |
| **2. Quality Control** | Assess read quality and length metrics | `NanoPlot`, `NanoStat`, `Seqkit` |
| **3. Alignment** | Map reads to reference genome | `Minimap2`, `Samtools` |
| **4. Variant Calling** | Detect SNVs and InDels | `Clair3` |
| **5. Annotation** | Predict and filter variant effects | `SnpEff`, `SnpSift` |  
## 2. Pipeline Tools

| Tool | Description | Resource Link
|------|--------------|--------|
| **NanoPlot** | Generates visual and statistical summaries of ONT read quality. | [NanoPlot](https://github.com/wdecoster/NanoPlot)
| **NanoStat** | Computes sequencing run metrics (read count, length, quality). | [NanoStat](https://github.com/wdecoster/nanostat)
| **Seqkit** | Fast toolkit for FASTA/FASTQ manipulation and statistics. | [Seqkit](https://github.com/shenwei356/seqkit)
| **Minimap2** | Rapid and accurate aligner for ONT long reads. | [Minimap2](https://github.com/lh3/minimap2)
| **Samtools** | Tools for SAM/BAM file manipulation and indexing. | [Samtools](https://github.com/samtools/samtools)
| **Clair3** | Deep learning–based variant caller optimized for ONT data. | [Clair3](https://github.com/nanoporetech/Clair3)
| **SnpEff** | Predicts functional consequences of variants using genomic annotations. | [SnpEff](https://pcingola.github.io/SnpEff/)
| **SnpSift** | Filters and annotates variants; integrates dbSNP/ClinVar data. | [SnpSift](https://pcingola.github.io/SnpEff/#snpsift)  
## 3. Installation Guide

Set up the environment using **WSL (Ubuntu)** and manage dependencies via **Conda** or **Python virtual environments**.  
Refer to each tool’s official documentation for detailed installation and version compatibility.

### A. WSL Setup
Enable **WSL** and install **Ubuntu** from Microsoft Store:
```bash
sudo apt update && sudo apt upgrade -y
sudo apt install default-jre -y  # Required for SnpEff & SnpSift
```

### B. Miniconda Setup

Download and Install Miniconda:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Add Conda Channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
### C. Tools Installation  
| Tool                    | Installation Method                                                                                                                                     |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Samtools**            | `conda install -c bioconda samtools`                                                                                                                    |
| **Seqkit**              | `conda install -c bioconda seqkit`                                                                                                                      |
| **Minimap2**            | `conda install -c bioconda minimap2`                                                                                                                    |
| **NanoPlot / NanoStat** | Create and activate virtual environment:<br>`python3 -m venv ~/nanoplot_env`<br>`source ~/nanoplot_env/bin/activate`<br>`pip install NanoPlot NanoStat` |
| **SnpEff / SnpSift**    | Download and configure from [SnpEff official site](https://pcingola.github.io/SnpEff/)                                                                  |
| **Clair3**              | Install via Bioconda: see [Clair3 Documentation](https://github.com/HKU-BAL/Clair3)                                                                     |  

### D. Reference Genome and Databases

Reference Genome (GRCh38):
Download from the NCBI Human Genome Resource. 
 ```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
minimap2 -d ref.mmi ref.fa # Index reference genome with minimap2
```
ClinVar Database (VCF):
```
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_latest.vcf.gz
```
### 4. Pipeline Commands

>Note: Threads are set to 10 for parallel processing — adjust according to CPU availability.

A. Merge Reads
```
cat *.fastq.gz > merged_reads.fastq.gz
```
B. Quality Control
```
source ~/nanoplot/bin/activate # Activate NanoPlot Environment
NanoPlot --fastq merged_reads.fastq.gz -o nanoplot_output --loglength
```
C. Alignment

Check Read Headers:
```
zcat merged_reads.fastq.gz | head -n 4
```
Align Reads Using Minimap2:
```
minimap2 -ax map-ont -t 10 -R '@RG\tID:FlowCell\tPL:ONT\tSM:Sample' ref.fa merged_reads.fastq.gz > aligned.sam
```
D. BAM Conversion, Sorting & Indexing
```
samtools view -@ 10 -Sb aligned.sam > aligned.bam
samtools sort -@ 10 -o sorted.bam aligned.bam
samtools index sorted.bam
```
E. Variant Calling (Clair3)
```
./run_clair3.sh \
--bam_fn sorted.bam \
--ref_fn ref.fa \
--threads 10 \
--platform ont \
--model_path /path/to/model \
--output clair3_output
```
F. Variant Annotation (SnpEff / SnpSift)

Add rsIDs and Annotate Variants:
```
# Clinvar rsIDs:
java -Xmx10g -jar SnpSift.jar annotate -v clinvar_latest.vcf.gz clair3_output.vcf > rsIDs.vcf
# Human Genome Database Annotation:
java -Xmx10g -jar snpEff.jar -v GRCh38_latest rsIDs.vcf > annotated.vcf
```
Filter Variants with HIGH or MODERATE Impact:
```
java -Xmx10g -jar SnpSift.jar filter \
"((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" annotated.vcf > filtered.vcf
```
### 5. References
- [NanoPlot Documentation](https://github.com/wdecoster/NanoPlot)  
- [NanoStat Documentation](https://github.com/wdecoster/nanostat)  
- [Seqkit Documentation](https://github.com/shenwei356/seqkit)  
- [Minimap2 Documentation](https://github.com/lh3/minimap2)  
- [Samtools Documentation](https://github.com/samtools/samtools)  
- [Clair3 Documentation](https://github.com/nanoporetech/Clair3)  
- [SnpEff & SnpSift Documentation](https://pcingola.github.io/SnpEff/)  
