---
title: Practical
---

## Finding a known genome in mixed microbial community

The following practical simulates the situation when we know what we are looking for and utilise a whole-genome alignment approach to find a particular bacterial strain or virus genome in a highly complex microbial community. This technique is especially useful if the known species (potential pathogen) is hard to detect / culture with traditional laboratory methods, amplification-based methods are not specific enough (e.g., often false positive due to a common species with highly similar genome), the laboratory method requires more time and/or more expensive.

::: {.callout-tip}
#### Key Points to consider

- The closer the genome you use to the known species / strain the higher sensitivity and specificity can be achieved. If you are monitoring an outbreak with this approach, the best results can be achieved if you can isolate the microbe at least once and perform a whole genome sequencing and de novo assembly on it.
- If you plan to quantify the tracked species, you can have a broad idea (relative abundance) by comparing the aligned read number to the total number of reads in the raw data. You can achieve much better quantification (both relative and absolute) if you spike in your sample (before DNA extraction) with a known bacteria or virus, using a well-defined amount. Ideally the spiked in species should be a distant species in terms of phylogeny and should have similarish genome size.
- If you don't have your own reference genome, try to find one in public databases that is potentially the closest to your geographical location but also a recent isolate.
:::

### Standard quality control and pre-processing of shotgun metagenomics raw data

Before we perform any analysis on the raw data it is important to perform the basic quality control checks and if needed certain pre-processing and filtering steps to ensure that we are working with high quality data. When you open a new terminal in the training environment, your working directory should be the `~/Course_Materials` folder. You can always check where you are in the filesystem by using the `pwd` command or just by checking your `bash` prompt.

Switch to the folder with the raw data, check the usage of the application and run `fastqc` on our first example raw data set (a pair of FASTQ files simulating a paired-end sequencing data).

```bash
cd sg_raw_data/
fastqc -h
fastqc mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz
```
::: {.callout-exercise}
#### Inspect the output of the FastQC application

FastQC generatesgraphical output report in `.html` format. This is often placed (and archived) together with the raw data files, so the quality measures can be quickly checked in the case of future usage.

Open the Html files and go through the graphs, discuss what you see. For a future reference, and to see more examples (good and bad data), please visit the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
:::

The next standard pre-processing step is to remove adapter, primer and other unwanted sequences from your reads. These sequence contents are side products of the next-generation sequencing technique, and as they are not coming from the template DNA, they can interfere from many downstream pipeline steps. One of the most commonly used tool for this purpose is [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Check the command line help for the application and run the filtering step on the raw sequencing data (please note, that you have to activate the `metagenomics` conda environment if it is not yet active).

```bash
conda activate metagenomics

cutadapt -h

cutadapt -a CTGTCTCTTATACACATCT -A ATGTGTATAAGAGACA \
-o mixed_bacterial_community_5M_noadapt_R1.fastq.gz -p mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz

```

The final standard pre-processing step is to remove bad quality sequencing from our raw data. Typically the 5' and 3' ends of the sequencing reads have lover quality (so it is worth to remove these ends), and it is also common that a few reads are in general bad quality. We will use the [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) application to perform this quality filtering step, please refer to the on-line website or the command line help (`trimmomatic -h`) for more information on usage and fine-tuning.

```bash
trimmomatic PE -phred33 mixed_bacterial_community_5M_noadapt_R1.fastq.gz mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixedcomm_forward_paired.fq.gz mixedcomm_forward_unpaired.fq.gz \
mixedcomm_reverse_paired.fq.gz mixedcomm_reverse_unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
::: {.callout-tip}
- Please note, that we are defining output file names for both the filtered paired-end and unpaired reads. While our input is coming from a paired-end library, it is possible that one mate of a read pair has bad quality while the other mate is good enough to carry forward. In this case the good read will be carried forward as unpaired read.
- To be able to use the `trimmomatic` tool efficiently, you have to have a good understanding of both the sequencing technique, the significance of potential errors and warnings listed in the `FastQC` report and the usage of the many different options and parameters in the `trimmomatic` application itself.
:::

::: {.callout-exercise}
#### Re-analyse the filtered raw data with FastQC and compare it to the original report
{{< level 2 >}}
Run the `fastqc` step again but this time giving the output files of the `trimmomatic` step to the application. Compare the reports and graphs and discuss the differences with the trainers.

::: {.callout-answer}
To run the `fastqc` on the `trimmomatic` output, you just simply list the new file names after the command. The new `.html` files will be generated next to the files:
```bash
fastqc mixedcomm_forward_paired.fq.gz \
mixedcomm_forward_unpaired.fq.gz \
mixedcomm_reverse_paired.fq.gz \
mixedcomm_reverse_unpaired.fq.gz
```
:::
:::

### Obtaing and preparing the reference genome for alignment

Reference genomes are commonly stored in [FASTA](21-pres.html#file-and-data-formats) formatted files. If the genome contains multiple DNA species (e.g., mutliple chromosomes or genome + plasmid) it is still stored in a single FASTA file, but in a multi-FASTA format, where each DNA species has its own FASTA header. These genome files can be generated (e.g., by bacterial colony sequencing and de novo assembly) or be downloaded from on-line genome databases. The most commonly used general databases are [Ensembl](https://www.ensembl.org/index.html), [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/), [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) and the metagenomics specific [MGnify](https://www.ebi.ac.uk/metagenomics).

There are multiple reference genome files saved in the `Course_Materials/sg_reference/` directory, for this exercise we will use the the reference genome of the Escherichia coli O157:H7 strain.

```bash
cd ~/Course_Materials/sg_reference/


```


```bash
####### Day 2, Alignment-based complete compositional profiling of complex microbial communities ########

# In our first example, we will use a custom build database to screen a complex microbial community
# Using traditional short read alignment and extracting coverage information

# As always we are starting the practical by being in the course_materials/ folder
cd sg_reference/

# Prepare the database and index it for biowtie2 usage
bowtie2-build -q mixed_bacterial_community_ncbi_genomes.fasta mixed_bacterial_community_ncbi_genomes.fasta
cd ..

# Perform the alignment of raw sequencing reads to the reference database set
bowtie2 -p 8 --fast --qc-filter -x sg_reference/mixed_bacterial_community_ncbi_genomes.fasta \
-1 sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz -2 sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz \
-S results/mixed_community_alignment.sam

# Let's look into the result .sam file and learn more about the .sam file format
less results/mixed_community_alignment.sam

# For long term storage, analysis, subsetting and visualisation (e.g., by BamView), it is better to store the alignment information in a binary compressed file format
# For this purpose, we use .bam file format.

# Let's convert our .sam file to a .bam file, and when doing it, sort it by position (most of the downstream analyses tools require sorted alignment files)
# (micromamba install -c bioconda samtools)
cd results/
samtools sort -@ 4 -O BAM -o mixed_community_alignment_sorted.bam mixed_community_alignment.sam

# Compare the file sizes of the original .sam and the binary and compressed .bam file
ls -lh

# Let's try the highest compression level to see how much we can reduce the file size
samtools sort -@ 4 -l 9 -O BAM -o mixed_community_alignment_sorted_smallest.bam mixed_community_alignment.sam
ls -lh

# Let's look at some statistics about the alignment
samtools flagstat mixed_community_alignment_sorted.bam

samtools stats mixed_community_alignment_sorted.bam > alignment_statistics.txt
less alignment_statistics.txt

samtools coverage mixed_community_alignment_sorted.bam

# Save this coverage file for future comparison
samtools coverage mixed_community_alignment_sorted.bam > mixed_community_coverage_from_alignment.txt
```

```bash
# Let's try a different method for compositional mapping from shotgun metagenomics data
metaphlan sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz,sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz \
--bowtie2out results/metaphlan.bowtie2.bz2 --nproc 5 --input_type fastq -o results/profiled_metagenome.txt

# Look into the results file
less results/profiled_metagenome.txt

# Keep only the full taxonomy lines for further analysis
grep "s__" results/profiled_metagenome.txt > results/profiled_metagenome_sp_only.txt

# R script TODO
```

```bash
####### Day 2, Direct analysis of short reads using K-mer and database-based methods ########

# Starting with mash, looking at the main functions
mash

# Use the mash info to find out details about sketch files
mash info sg_reference/RefSeqSketches.msh | head -n 20
mash info sg_reference/RefSeqSketchesDefaults.msh | head -n 20

# Check the help for sketching
mash sketch -h

# List the reference sequence folder and change to the directory
ls sg_reference/
cd sg_reference/

# We will create a new sketch file from mixed_bacterial_community_ncbi_genomes.fasta file
# First, let's see what is inside the file

less mixed_bacterial_community_ncbi_genomes.fasta
grep ">" mixed_bacterial_community_ncbi_genomes.fasta
grep ">" mixed_bacterial_community_ncbi_genomes.fasta | wc -l

# Create a "quick and light" sketch
mash sketch -i -s 500 -k 16 -o mixed_community mixed_bacterial_community_ncbi_genomes.fasta

# Compare the original fasta file size with the sketch size
ls -lh

# Check the contents of the sketch file
mash info mixed_community.msh

# Create a more specific and high kmer count number sketch from the same fasta file
mash sketch -i -s 5000 -k 21 -o mixed_community_hr mixed_bacterial_community_ncbi_genomes.fasta

# Compare the original fasta file size with the sketch size
ls -lh

# Check the contents of the sketch file
mash info mixed_community_hr.msh

# Check the running options and parameters for mash screen
mash screen -h

# Inspect the two NZ_ fasta files in the reference directory
less NZ_CP034931.fa
grep ">" NZ_CP034931.fa

less NZ_CP038419_1.fasta
grep -c ">" NZ_CP038419_1.fasta

# Screen the two fasta files in the reference directory
mash screen mixed_community.msh NZ_CP034931.fa
mash screen mixed_community.msh NZ_CP038419_1.fasta

# Run both again with the "Winner-takes-all" option and compare the results with the provious ones
mash screen -w mixed_community.msh NZ_CP034931.fa
mash screen -w mixed_community.msh NZ_CP038419_1.fasta

# Compare the results of screen using the smaller and the larger sketch set we created, with and without the -w option
mash screen mixed_community.msh NZ_CP038419_1.fasta
mash screen mixed_community_hr.msh NZ_CP038419_1.fasta

mash screen -w mixed_community.msh NZ_CP038419_1.fasta
mash screen -w mixed_community_hr.msh NZ_CP038419_1.fasta

# Screen our single fasta file with the NCBI reference collections
mash screen -w RefSeqSketches.msh NZ_CP034931.fa # DON'T ASK TO DO JUST DEMONSTRATE!!!
mash screen -w RefSeqSketches.msh NZ_CP034931.fa | sort -gr -k 1 | head
mash screen -w RefSeqSketchesDefaults.msh NZ_CP034931.fa | sort -gr -k 1 | head

# Screen the simulated shotgun sequencing data using both the high resoluton custom sketch file and the NCBI file
mash screen -w mixed_community_hr.msh ../sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz ../sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz
mash screen -w RefSeqSketchesDefaults.msh ../sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz ../sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz | sort -gr -k 1 | head -n 50

# Try to accelerate the process by using multiple CPU cores
time mash screen -w RefSeqSketchesDefaults.msh ../sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz ../sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz | sort -gr -k 1 | head -n 50
time mash screen -w -p 4 RefSeqSketchesDefaults.msh ../sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz ../sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz | sort -gr -k 1 | head -n 50

# To be able to run the next part, you need to install SRST2 in a new empty environment
# please replace "conda" command with the package manager you have (mamba / micromamba)
conda deacticvate
conda create -n srst2
conda activate srst2
conda install -c bioconda -c conda-forge srst2


# Let's see how srst2 works
srst2 -h

# Try to run SRST2 on our mixed community sequencing data
srst2 --input_pe sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz --output results/amr_genes.txt

# SRST2 will drop an error as the raw data naming is strict when using the algorithm, 
# forward sequencing data file name has to end as "_1.fastq.gz" and reverse data file name has to end as "_2.fastq.gz"

# To be sure that our filenames remain the original (to be compatible with previous commands)
# and to not duplicate data (and fill upt our storage quicker), we create symbolic links
# These are virtual files (with any given filename) that are pointing to an existing file
cd sg_raw_data/
ln -s mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_1.fastq.gz
ln -s mixed_bacterial_community_5M_R2.fastq.gz mixed_bacterial_community_5M_2.fastq.gz
cd ..

# Download the CARD database from the SRST2 github website to the sg_reference directory
# https://github.com/katholt/srst2/blob/master/data/CARD_v3.0.8_SRST2.fasta
# And index the fasta file with bowtie2
bowtie2-build CARD_v3.0.8_SRST2.fasta CARD_v3.0.8_SRST2.fasta

srst2 --input_pe sg_raw_data/mixed_bacterial_community_5M_1.fastq.gz sg_raw_data/mixed_bacterial_community_5M_2.fastq.gz --output results/amr_genes.txt --gene_db sg_reference/CARD_v3.0.8_SRST2.fasta



```