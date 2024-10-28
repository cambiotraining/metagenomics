---
title: De novo and Hi-C metagenomics (theory)
---
## De novo metagenomics and Hi-C metagenomics

The slides for this presentation section can be found on google slides on the [following link](https://docs.google.com/presentation/d/1IUciE3_K6XK1Xy6koYxMDPwr5kILl1GYRAaT0t-qPss/edit?usp=sharing)

### Fundamentals of de novo and Hi-C metagenomics

In the de novo approach of shotgun metagenomics data analysis, we treat the sample as the mixture of completely unknown genomes and use computational and mathematical algorithms to combine raw reads and build contigs (continuous DNA sequences) and then cluster contigs to form assemblies (representing potential individual genomes). Compared to the reference genome based alignment methods (where the composition of our reference genome / gene / kmer set is limiting what we can find), the de novo approach is mostly unbiased, in many cases we have the same chance to assemble a completely unknown commensal and a well-known E. coli.

::: {.callout-tip}
#### Disadvanteges compared to alignment based methods

- Sensitivity is much lower, low abundance species are much harder to detect and assemble
- Significantly higher computational resource need, several steps require high amount of CPU cores (>24) and RAM (>100GB). To analyse a few samples with the de novo approach, usually an in-house server is needed, for projects involving large number of samples usually access to high-performance computing (HPC) services is necessary.
- Genetically closely related species are usually hard to separate
- Besides the higher computational resource requirement, this approach requires significantly more time to run the bioinformatics steps and also more advanced bioinformatics knowledge
:::

### Quality control and pre-processing during de novo analysis

As the de novo approach uses the same shotgun sequencing raw data as the alignment based approach, most of the QC and pre-pocessing steps are the same as we discussed in the previous chapter. We use `FastQC` to check the raw data quality, `trimmomatic` and `cutadapt` to trim the raw reads and remove sequencing method introduced extra DNA content. There are two extra pre-processing steps usually performed before we start the next stage and these are both helping / enhancing the main assembly step.First, as the assembly step resource and time need is going up significantly with the amount of input data, it is recommended to remove any duplicates (PCR or optical duplicates) from the raw read data. This is usually performed by the `clumpify.sh` script from the `bbmap` software package. Second, the de novo assembly efficiency can be significantly increased (in terms of e.g., reducing the fragmentation of contigs) by introducing longer sequencing reads. Short read sequencing techniques have their length limitations, but we can still make "longer" reads if we merge those read pairs that have overlapping ends (if the sequenced DNA fragment is shorter than the length of the sequencing of the read pair). This step is usually done by the application `FLASh` obtainable through [conda](https://anaconda.org/conda-forge/flash) or from the [original website](https://ccb.jhu.edu/software/FLASH/).

### De novo assembly step

The pre-processed and quality filtered data is fed into the assembly process. During the de novo assembly the software application fragments the short reads (usually in multiple round into different sized k-mers) and try to find matching fragments. The matching fragments then will be used to find reads that are partially overlapping, so most probably originated from the same template (same genome). The aim of the assembly step is to build up as long continuous DNA fragments (contigs) as possible by combining overlapping short reads. Ideally a whole genome should be built up as a single contig (by reaching the starting point of the contig in circular genomes) but this is very rarely the case. Sequence elements that are found in multiple copies in a single genome or between multiple genomes (intra and inter-genomics repeats) usually break the contig elongation as these repeats have multiple "exit" sequence pathways. The two most popular metagenomics de novo assemblers are [metaSPAdes](https://github.com/ablab/spades) and [MEGAHIT](https://github.com/voutcn/megahit).

### Binning step

The result of the de novo assembly step is one big multi-FASTA file containing 10s or often 100s of thousands variable sized contigs. These are all shorter or longer fragments of a genome or a plasmid, but at this point everything is in the same cluster in one single file. The aim of the binning step is to cluster contigs to groups (bins) that are potentially representing individual genomes. This step is very challenging and yet the most unreliable part of the de novo metagenomics pipeline. There are three different approaches we can perform binning with all have their pros and cons.

Supervised binning, when we use reference genomes to find matching contigs:

- <span style="color:green">**(+)**</span> Fast, reliable and simple
- <span style="color:green">**(+)**</span> Works well for known genomes
- <span style="color:red">**(-)**</span> High dependency on prior knowledge
- <span style="color:red">**(-)**</span> Knowledge bias can introduce composition bias
- <span style="color:red">**(-)**</span> Works poorly for less known microbiomes

Unsupervised binning, when we use genomics properties (e.g., GC content, coverage) to find matching contigs:

- <span style="color:green">**(+)**</span> No prior knowledge needed
- <span style="color:green">**(+)**</span> Good for less studied microbiomes
- <span style="color:red">**(-)**</span> Less sensitivity and specificity for known microbiomes
- <span style="color:red">**(-)**</span> Usually slow and complex procedure
- <span style="color:red">**(-)**</span> Problems separating closely related species
- Works better with multiple samples in parallel

Using long read sequencing data to find contigs from the same genome:

- <span style="color:green">**(+)**</span> No prior knowledge needed
- <span style="color:green">**(+)**</span> Can detect segmental rearrangements
- <span style="color:red">**(-)**</span> Expensive and requires special laboratory procedures
- <span style="color:red">**(-)**</span> Lack of hybrid (short and long read) pipelines and applications
- <span style="color:red">**(-)**</span> Method specific issues with accuracy

### Hi-C metagenomics

Performing Hi-C metagenomics requires an extra short-read sequencing library besides the standard shotgun metagenomics sequencing. The extra library includes special steps during the laboratory protocol and extra computational steps in the analysi pipeline. The main advantage of the Hi-C metagenomics approach that we create physical links within the cell (cross-links using formaldehyde before cell lysis) between random parts of the genome and between genome and extra-chromosomal DNA (e.g., plasmids). The cross-link information created by the laboratory protocol than will be used during the computational pipeline to help binning (contigs will have physical evidence of belonging to the same cell) and in associating extra chromosomal DNA to the host cell. Identification of the host species for plasmids is a unique feature of Hi-C metagenomics and cannot be achieved by any other metagenomics technique. This makes the Hi-C metagenomics an extremely useful tool in the research / diagnostics / monitoring of antimicrobial resistance presence, transmission and spreading.