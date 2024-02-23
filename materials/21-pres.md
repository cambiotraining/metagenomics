---
title: Presentations
---

## Data formats, QC and data management

The slides for this presentation section can be found on google slides on the [following link](https://docs.google.com/presentation/d/1q8iTHEGhUeIv_D8RKj-NycbKrC-YsPvd-Lup4oDi-kU/edit?usp=sharing)

### File and data formats

During the analysis of shotgun metagenomics data you will work with different data types. Most of these are plain text files (usually compressed by `gzip`), but having special syntax and internal structure. Understanding how these files look like and especially what purpose they are serving is an important basic knowledge. The most commonly used file formats are listed below, their structure and syntax are detailed in the presentation material.

| File type  | File extension                 | Main purpose                                                                                                                                                                       |
|------------|--------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FASTA | .fa, .fasta, .fna              | Storing nucleotide or protein sequence with the sequence identifier (e.g., gene name)                                                                                              |
| FASTQ | .fastq, .fq, .fastq.gz, .fq.gz | Storing sequencing raw data, the sequence is stored together with the acquisition quality                                                                                          |
| SAM   | .sam                           | Storing alignment data, keeping both the sequence and quality information from the FASTQ file but also storing information on the genomic position that the sequence is aligned to |
| BAM   | .bam                           | The binary compressed format of the SAM file, saves significant amount of storage space compared to SAM file                                                                       |
| VCF   | .vcf                           | Storing genetic variance information based on the position of the reference genome                                                                                                 |

::: {.callout-tip}
#### Key Points

- While there are standard tools to handle the above mentioned file types / formats, so most of the time you don't need to open and inspect them, for troubleshooting purposes it is always good to know more about the format and syntax of these files.
- If you convert a FASTQ file to a FASTA file (as your application can only read in FASTA format), be prepared that you will loose the sequencing quality information. If you need to do this, it is recommended to perform a quality based filtering first and only convert those sequence records that passed the filter.
- You can save significant amount of storage space by compressing (by `zip` or `gzip`) FASTA, FASTQ and VCF files.
- You can save significant amount of storage spae if you convert your SAM files to BAM files (they store exactly the same information)
:::

### Quality control of raw data

Every bioinformatics pipeline that processes next-generation sequencing data has to start with quality assessment. There are several measurable properties that users can extract from the raw data (FASTQ files), and these measures can indicate systematic errors that were introduced during sample collection, adapter ligation, library creation or even during the sequencing run. It is important to check all the raw data files, it may turn out that certain problems are only coming up with certain samples or certain batches. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is the most commonly used application for raw data quality assessment. 

::: {.callout-tip}
#### Tip

It is a good practice to archive raw data files with the permanent output file of the FastQC tool (usually a `.html` report file). This way the quality of the raw data is always available even years later when the data is re-used.
:::

### Quality control of the bioinformatics pipeline

Bioinformatics pipelines are diverse and the observed changes, transformations of the data can have different meaning even for the same application in different circumstances. E.g., finding high amount of host DNA contamination during a metagenomics pipeline in faecal samples may raise concerns about the DNA extraction protocol (potential bacterial cell lysis problems), while the similar high amount of host DNA in a metagenomics sample that is coming from swabbing (potentially scraping the host epithelium more) can be normal. Understanding the pipeline steps together with being aware of the potential experimental biases is a very important factor in controlling the quality of a metagenomics pipeline. There are tools to help extracting meaningful information from log files (e.g., [MultiQC](https://multiqc.info)), but always adapt the "normal ranges" to your own experimental setup.

### Data management

Managing data from metagenomics experiments can be a challenge due to the size of raw data, intermediate and processed files. During the planning phase of large metagenimics projects (especially long-term metagenomics based services), it is crucial to calculate with the expenses on data storage, archiving solutions and also on computational capacity.

::: {.callout-tip}
#### Key Points

- __Data in one copy is not data__, always keep a copy of your most valuable files (usually the raw data)
- Don't rely on public repositories as data storage or backup solution. These providers don't take any responsibility for your data, it serves as a platform for sharing your results with the public, but not as an archiving solution.
- Be aware of your country's policies and legislation on handling human samples. Even a human faecal sample contains significant amount of host epithelial cells and as a consequence human DNA will be sequenced. Consider "de-humanise" your raw data (remove short reads that are aligning to human genome) in an early step of your pipeline. Depositing metagenomics sample with human DNA information in it is also more complicated.
- Delete intermediate files as soon as you can, keep only those that required a lot of computational resources (so you don't have to re-do them).
- Compress files to save space, most of the text based formats can be efficiently compressed by the common `zip` or `gzip` algorithms. Alignment data should always be stored in BAM files instead of the SAM format.

:::

[Alignment based metagenomics](https://docs.google.com/presentation/d/1Q1mzNO-nOi0dYjK1l9Ig6b1Gpu-JdRvuQkxvXe-ozJU/edit?usp=sharing)