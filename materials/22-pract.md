---
title: Practical
---

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