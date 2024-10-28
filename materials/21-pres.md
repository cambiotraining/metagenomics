---
title: QC and alignment based metagenomics (theory)
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

## Alignment based metagenomics

The slides for this presentation section can be found on google slides on the [following link](https://docs.google.com/presentation/d/1Q1mzNO-nOi0dYjK1l9Ig6b1Gpu-JdRvuQkxvXe-ozJU/edit?usp=sharing)

### Fundamentals of alignment based approaches

In the alignment based metagenomics approaches we try to find the best matches for our sequencing data among known reference genomes. Match can mean perfect (nucleotide by nuclotide), imperfect (nuclotide changes, insertions, deletions allowed) alignments along the full length of the short read (global alignment) or just within fragments of it (local alignment, hash based methods). From the alignment results, we can derive multiple information: (i) presence / absence of a certain species; (ii) relative abundance of a species (based on alignment coverage data); (iii) genetic distance of a certain species / strain from the closest homologue (based on genetic variations and genome structural homology); (iv) presence / absence of certain genes in a species (e.g., antimicrobial resistance or virulance genes).

::: {.callout-tip}
#### Key Points

- You can only find species / genes with alignment based techniques if you have something similar in your reference database.
- Due to our significantly better knowledge on pathogen species, alignement based methods will always be a bit biased towards well-known species (compared to e.g., commensals in the gut microbiome).
- Alignment based methods are very sensitive and specific to find known genetic material but have severe limitations in discovering emerging, previously unknown species / plasmids / genes.

:::

### K-mer based profiling

In k-mer based methods we randomly extract short nucleotide sequences from reference genomes (fragment size is usually between 16-32 nucleotides, number of fragments ranging from few hundred to several thousands) and screen raw sequencing data (FASTQ files) or assembled contigs (FASTA files) for perfect matching. The theory is, that if we have a certain genome / species in our mixed community high percentage of the fragments, originally extracted from that genome, will show perfect matching with the short reads or the assembled contigs. By defining the sensitivity threshold, we can obtain presence / absence information, by counting an average number of exact matches for all k-mers from the same genome, we can derive abundance information. The most commonly used general k-mer based application is [Mash](https://mash.readthedocs.io/en/latest/), while a bit more microbiology and metagenomics optimised application is [Kraken](https://ccb.jhu.edu/software/kraken2/).

::: {.callout-tip}
#### Key Points

- K-mer based methods provide the fastest way to profile shutgun metagenomics raw data. These methods have also relatively low computational resource requirements (CPU, memory).
- As all alignment based methods, k-mer based methods sensitivity and reliability is highly dependent on our previous knowledge (our initial database).
- While species level identification is possible with k-mer based approach, strain specific resolution is challenging.
- As the selection of k-mers from the reference genomes is random, the presence / absence of specific k-mers rarely have any biological meaning.
:::

### Marker gene detection based profiling

As our knowledge grows in identifying novel microbial genomes, we can compare these genomes and identify genes / gene families that are specific to certain taxonomic levels. This way, we are able to use small subsets of the genome (smaller database, faster bioinformatics pipeline) in the identification of our complex community members. The method provides a biologically meaningful way (as we actually searching for genes) to detect the presence of a previously completely unknown genome and potentially assign to a taxonomic category (even if not identified on species level, we may be able to assign a genus or a family). In this approach the reference gene databases are provided by the software developers (as their generation / modification is not a simple task), so the user has to be very careful of using a well maintained and up-to-date application / database.

::: {.callout-tip}
#### Key Points

- Marker gene based methods are usually require slightly more computational resources compared to the k-mer based methods but still have relatively low resource needs (a high-spec desktop or laptop can usually handle the task).
- Marker gene based methods can reliably detect the microbial community members on the species level, further refinement of the gene sequences and polymorphisms (e.g., by using [StrainPhlAn](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1)) can resolve the genomes in strain level.

:::

### Whole genome alignment based detection and profiling

The most sensitive and specific detection of a certain species or strain of a microbial genome from shotgun metagenomics sequencing can be performed by aligning all raw reads to known whole genomes. While it is almost impossible to use all the known microbial genomes in one step and align reads to them, after a k-mer or marker gene based discovery, the high confidence / abundance genomes can be extracted from the raw data by specifically aligning the reads to the known genomes. Alignment based methods can also be utilised to reduce the amount of raw sequencing reads (for e.g., de novo metagenomics assembly where the amount of reads is in strong positive correlation with the analysis time and memory requirement), highly abundant known genomes and/or the host genome can be used to remove reads that are aligning to these and use the "remaining" raw data to do de novo assembly and novel genome discovery. An extreme example for the usage of alignment based methods is to eliminate certain genomes from a low complexity community (e.g., human blood where we suspect an unknown pathogen), and enriching the raw data for reads originating from the unknown genome.

::: {.callout-tip}
#### Key Points

- Whole genome alignment methods in metagenomics rarely used to map the full composition landscape of the community (exceptions are the artificial communities with known components).
- Use full genome alignment approach if you would like to find a specific species / strain with the highest specificity and sensitivity.
- Use the alignment based in combination with a de novo approach to to save computational time and be able to use less resources in the de novo assembly.
- Utilise the whole genome alignment approach to eliminate unwanted and/or well-known genomes from your raw data to make the furhter analysis steps more efficient.

:::