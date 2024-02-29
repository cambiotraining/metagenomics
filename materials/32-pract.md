---
title: Practical
---

```bash
########  Day 2 practical 2: Finding an unknown genome in a low complexity microbial community #######

# In this practical, we simulate a metagenomics approach when we try to find and assemble a completely unknown genome
# We sequence a low diversity community (in this case we only have "unwanted" human sequncing background)

# Be sure you are in the "metagenomics" environment
# You may need
conda deactivate
conda activate metagenomics

# Let's look at the help of one of the most commonly used aligner algorithm
bowtie2 -h

# We need to index our database first
cd sg_reference/
bowtie2-build -q Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa

# Notice the new files that were created during indexing
ls
cd ..

# Align the raw reads to the Human chromosome 22 reference sequence
bowtie2 --qc-filter -p 8 --local -x sg_reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa -1 sg_raw_data/unknown_pathogen_R1.fastq -2 sg_raw_data/unknown_pathogen_R2.fastq -S results/unknown_pathogen.sam

# Inspect the result .sam file
less results/unknown_pathogen.sam

# For this analysis type we don't really need the alignment result to the human genome, we know it has to contain human DNA, but we want to filter it out to enrich for reads from the unknown pathogen
# Let's re-run bowtie2 with some modifications
bowtie2 --qc-filter -p 8 --local -x sg_reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa \
-U sg_raw_data/unknown_pathogen_R1.fastq,sg_raw_data/unknown_pathogen_R2.fastq --un sg_raw_data/unknown_enriched.fastq > /dev/null

# We needed to input the raw reads as unpaired as there were some error during the read-sett simulation and it doesn't align well to the reference

# In the final step, as we expect the remaining DNA is coming from a single source, we use traditional de novo assembly to reconstruct the genome
# We use the same software / method to build bacterial genomes from colony sequencing

# Check how to use the spades assembler
spades.py

# Run the assembly using all the default settings, only giving the input raw data file and the output folder
spades.py -s sg_raw_data/unknown_enriched.fastq -o results/unknown_genome

# The assembly results will be in the results/unknown_genome/ folder, check the output files, logs, warnings
cd results/unknown_genome/

```

```bash
########  Day 3 practical 1: De novo metagenomics practical  #######

# The following pre-processing needs to be done on short read sequencing data for both de novo and alignment based methods

cd sg_raw_data/

# Running FastQC for checking the quality of sequencing
fastqc -h

fastqc mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz

# If the error "java: symbol lookup error: java: undefined symbol: JLI_StringDup" occurs, do the following

micromamba deactivate
micromamba create -n fastqc
micromamba activate fastqc
micromamba install -c bioconda fastqc

fastqc mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz

# New .zip and .htm files were created with the raw data file names, inspect the created html file in a browser
 
# Dea-activate the fastqc environment and go back to the metagenomics env
micromamba deactivate
micromamba activate metagenomics

# Even if we didn't see significant amount of adapter sequences on the FastQC report, still it is a standard step to remove these (as usually they exist in higher abundance)
cutadapt -a CTGTCTCTTATACACATCT -A ATGTGTATAAGAGACA \
-o mixed_bacterial_community_5M_noadapt_R1.fastq.gz -p mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz

# Next step is to remove bad quality reads and nucleotides from our raw read sequences
trimmomatic -h

trimmomatic PE -phred33 mixed_bacterial_community_5M_noadapt_R1.fastq.gz mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixedcomm_forward_paired.fq.gz mixedcomm_forward_unpaired.fq.gz \
mixedcomm_reverse_paired.fq.gz mixedcomm_reverse_unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# again if you get the above java error, do the following steps:
micromamba deactivate
micromamba create -n trimmomatic
micromamba activate trimmomatic
micromamba install -c bioconda trimmomatic

# And run the trimmomatic step again

# Don't forget to go back to the metagenomics env after you finish
micromamba deactivate
micromamba activate metagenomics

# In the next step we will merge overlapping reads (like we did in the 16S pipeline)
# The main aim of this is to reduce the amount of input sequence and to create longer merged reads that can help during the de novo assembly

flash --help

flash mixedcomm_forward_paired.fq.gz mixedcomm_reverse_paired.fq.gz
ls -ltr

# We will carry forward 3 files for the de novo assembly, two for the paired end reads and one for unpaired and merged reads
# Let's merge in the unpaired reads (came from the trimmomatic step) to the merged reads file
# Please be aware that the merged file is plain fastq file, while the unpaired files are gzip-ed, so needed to be unzipped

zcat mixedcomm_forward_unpaired.fq.gz mixedcomm_reverse_unpaired.fq.gz >> out.extendedFrags.fastq

# Perform the de novo assembly step
metaspades.py

metaspades.py -t 8 -1 out.notCombined_1.fastq -2 out.notCombined_2.fastq -s out.extendedFrags.fastq -o ../results/mixed_comm_de_novo


```
